# -*- coding: utf-8 -*-
"""

.. codeauthor:: Rasmus Handberg <rasmush@phys.au.dk>
"""

import logging
import sys
import os.path
import glob
from astropy.table import Table
import hashlib
import shutil
import gzip
from zipfile import ZipFile
import tempfile
import re
from flows.aadc_db import AADC_DB
from flows.plots import plt, plot_image
from flows.load_image import load_image

#--------------------------------------------------------------------------------------------------
def flows_get_archive_from_path(fname, archives_list=None):
	"""
	Translate full path into AADC archive identifier and relative path.

	It is highly recommended to provide the list with that call
	to this function since it will involve a query to the database
	at every call.
	"""

	archive = None
	relpath = None

	# Get list of archives, if not provided with call:
	if archives_list is None:
		with AADC_DB() as db:
			db.cursor.execute("SELECT archive,path FROM aadc.files_archives;")
			archives_list = db.cursor.fetchall()

	# Make sure folder is absolute path
	folder = os.path.abspath(fname)

	# Loop through the defined archives and find one that matches:
	for opt in archives_list:
		archive_path = opt['path']
		if archive_path is not None and archive_path != '':
			archive_path = archive_path.rstrip('/\\') + os.path.sep
			if folder.startswith(archive_path):
				archive = int(opt['archive'])
				relpath = folder[len(archive_path):].replace('\\', '/')
				break

	# We did not find anything:
	if archive is None:
		raise Exception("File not in registred archive")

	return archive, relpath

#--------------------------------------------------------------------------------------------------
def get_filehash(fname):
	"""Calculate SHA1-hash of file."""
	buf = 65536
	s = hashlib.sha1()
	with open(fname, 'rb') as fid:
		while True:
			data = fid.read(buf)
			if not data:
				break
			s.update(data)

	sha1sum = s.hexdigest().lower()
	if (len(sha1sum) != 40): raise Exception("Invalid file hash")
	return sha1sum

#--------------------------------------------------------------------------------------------------
def optipng(fpath):
	os.system('optipng -preserve "%s"' % fpath)

#--------------------------------------------------------------------------------------------------
def create_plot(filepath):

	output_fpath = os.path.abspath(re.sub(r'\.fits(\.gz)?$', '', filepath) + '.png')

	img = load_image(filepath)

	cmap = plt.cm.Blues
	cmap.set_bad('w')

	fig = plt.figure(figsize=(12,12))
	ax = fig.add_subplot(111)
	plot_image(img.clean, ax=ax, scale='linear', percentile=[5, 99], cmap=cmap, cbar='right')
	fig.savefig(output_fpath, bbox_inches='tight')
	plt.close(fig)

	optipng(output_fpath)

#--------------------------------------------------------------------------------------------------
def ingest_from_inbox():

	rootdir_inbox = '/flows/inbox'
	rootdir = '/flows/archive'

	logger = logging.getLogger(__name__)

	# Check that root directories are available:
	if not os.path.isdir(rootdir_inbox):
		raise FileNotFoundError("INBOX does not exists")
	if not os.path.isdir(rootdir):
		raise FileNotFoundError("ARCHIVE does not exists")

	with AADC_DB() as db:
		# Get list of archives:
		db.cursor.execute("SELECT archive,path FROM aadc.files_archives;")
		archives_list = db.cursor.fetchall()

		for inputtype in ('science', 'templates', 'subtracted'):
			for fpath in glob.iglob(os.path.join(rootdir_inbox, '*', inputtype, '*')):
				logger.info("="*72)
				logger.info(fpath)

				# Find the uploadlog corresponding to this file:
				db.cursor.execute("SELECT logid FROM flows.uploadlog WHERE uploadpath=%s;", [os.path.relpath(fpath, rootdir_inbox)])
				row = db.cursor.fetchone()
				if row is not None:
					uploadlogid = row['logid']
				else:
					uploadlogid = None

				# Only accept FITS file, or already compressed FITS files:
				if not fpath.endswith('.fits') and not fpath.endswith('.fits.gz'):
					if uploadlogid:
						db.cursor.execute("UPDATE flows.uploadlog SET status='Invalid file type' WHERE logid=%s;", [uploadlogid])
						db.conn.commit()
					logger.error("Invalid file type: %s" % os.path.relpath(fpath, rootdir_inbox))
					continue

				# Get the name of the directory:
				# Not pretty, but it works...
				target_dirname = fpath[len(rootdir_inbox)+1:]
				target_dirname = target_dirname.split(os.path.sep)[0]

				# Convert directory name to target
				db.cursor.execute("SELECT targetid,target_name FROM flows.targets WHERE target_name=%s;", [target_dirname])
				row = db.cursor.fetchone()
				if row is None:
					logger.error('Could not find target: %s', target_dirname)
					continue
				targetid = row['targetid']
				targetname = row['target_name']

				if not fpath.endswith('.gz'):
					# Gzip the FITS file:
					with open(fpath, 'rb') as f_in:
						with gzip.open(fpath + '.gz', 'wb') as f_out:
							f_out.writelines(f_in)

					# We should now have a Gzip file instead:
					if os.path.isfile(fpath) and os.path.isfile(fpath + '.gz'):
						# Update the log of this file:
						if uploadlogid:
							db.cursor.execute("UPDATE flows.uploadlog SET uploadpath=%s WHERE logid=%s;", [os.path.relpath(fpath+'.gz', rootdir_inbox), uploadlogid])
							db.conn.commit()

						os.remove(fpath)
						fpath += '.gz'
					else:
						raise Exception("Gzip file was not created correctly")

				if inputtype == 'science':
					newpath = os.path.join(rootdir, targetname, os.path.basename(fpath))
					datatype = 1
				elif inputtype == 'templates':
					newpath = os.path.join(rootdir, targetname, inputtype, os.path.basename(fpath))
					datatype = 3
				elif inputtype == 'subtracted':
					newpath = os.path.join(rootdir, targetname, inputtype, os.path.basename(fpath))
					datatype = 4

					original_fname = os.path.basename(fpath).replace('diff.fits', '.fits')
					db.cursor.execute("SELECT fileid FROM flows.files WHERE targetid=%s AND datatype=1 AND path LIKE %s;", [targetid, '%/' + original_fname])
					subtracted_original_fileid = db.cursor.fetchone()
					if subtracted_original_fileid is None:
						if uploadlogid:
							db.cursor.execute("UPDATE flows.uploadlog SET status='original science image not found' WHERE logid=%s;", [uploadlogid])
							db.conn.commit()
						logger.error("ORIGINAL SCIENCE IMAGE COULD NOT BE FOUND: %s" % os.path.basename(fpath))
						continue
					else:
						subtracted_original_fileid = subtracted_original_fileid[0]
				else:
					raise Exception("Not understood, Captain")

				logger.info(newpath)

				if os.path.exists(newpath):
					logger.error("Already exists")
					if uploadlogid:
						db.cursor.execute("UPDATE flows.uploadlog SET status='Already exists: file name' WHERE logid=%s;", [uploadlogid])
						db.conn.commit()
					continue

				archive, relpath = flows_get_archive_from_path(newpath, archives_list)

				db.cursor.execute("SELECT fileid FROM flows.files WHERE archive=%s AND path=%s;", [archive, relpath])
				if db.cursor.fetchone() is not None:
					print("ALREADY DONE")
					continue

				# Calculate filehash of the file being stored:
				filehash = get_filehash(fpath)

				# Check that the file does not already exist:
				db.cursor.execute("SELECT fileid FROM flows.files WHERE filehash=%s;", [filehash])
				if db.cursor.fetchone() is not None:
					print("ALREADY DONE: Filehash")
					if uploadlogid:
						db.cursor.execute("UPDATE flows.uploadlog SET status='Already exists: filehash' WHERE logid=%s;", [uploadlogid])
						db.conn.commit()
					continue

				try:
					img = load_image(fpath)
				except:
					logger.exception("Could not load FITS image")
					continue

				if img.site['siteid'] is None:
					logger.error("Unknown SITE")
					continue

				try:
					# Copy the file to its new home:
					os.makedirs(os.path.dirname(newpath), exist_ok=True)
					shutil.copy(fpath, newpath)

					filesize = os.path.getsize(fpath)

					if not fpath.endswith('-e00.fits'):
						create_plot(newpath)

					db.cursor.execute("INSERT INTO flows.files (archive,path,targetid,datatype,site,filesize,filehash,obstime,photfilter,exptime,available) VALUES (%(archive)s,%(relpath)s,%(targetid)s,%(datatype)s,%(site)s,%(filesize)s,%(filehash)s,%(obstime)s,%(photfilter)s,%(exptime)s,1) RETURNING fileid;", {
						'archive': archive,
						'relpath': relpath,
						'targetid': targetid,
						'datatype': datatype,
						'site': img.site['siteid'],
						'filesize': filesize,
						'filehash': filehash,
						'obstime': img.obstime.mjd,
						'photfilter': img.photfilter,
						'exptime': img.exptime
					})
					fileid = db.cursor.fetchone()[0]

					if inputtype == 'subtracted':
						db.cursor.execute("INSERT INTO flows.files_cross_assoc (fileid,associd) VALUES (%s,%s);", [fileid, subtracted_original_fileid])

					if uploadlogid:
						db.cursor.execute("UPDATE flows.uploadlog SET fileid=%s,status='ok' WHERE logid=%s;", [fileid, uploadlogid])

					db.conn.commit()
				except:
					db.conn.rollback()
					if os.path.exists(newpath):
						os.remove(newpath)
					raise
				else:
					logger.info("DELETE THE ORIGINAL FILE")
					if os.path.isfile(newpath):
						os.remove(fpath)
					if uploadlogid:
						db.cursor.execute("UPDATE flows.uploadlog SET uploadpath=NULL WHERE logid=%s;", [uploadlogid])
						db.conn.commit()

#--------------------------------------------------------------------------------------------------
def ingest_photometry_from_inbox():

	rootdir_inbox = '/flows/inbox'
	rootdir_archive = '/flows/photometry'

	logger = logging.getLogger(__name__)

	# Check that root directories are available:
	if not os.path.isdir(rootdir_inbox):
		raise FileNotFoundError("INBOX does not exists")
	if not os.path.isdir(rootdir_archive):
		raise FileNotFoundError("ARCHIVE does not exists")

	with AADC_DB() as db:
		# Get list of archives:
		db.cursor.execute("SELECT archive,path FROM aadc.files_archives;")
		archives_list = db.cursor.fetchall()

		# Delete left-over files in the tables, that have been removed from disk:
		db.cursor.execute("SELECT logid,uploadpath FROM flows.uploadlog WHERE uploadpath IS NOT NULL;")
		for row in db.cursor.fetchall():
			if not os.path.isfile(os.path.join(rootdir_inbox, row['uploadpath'])):
				print("DELETE THIS FILE: " + row['uploadpath'])
				db.cursor.execute("DELETE FROM flows.uploadlog WHERE logid=%s;", [row['logid']])
				db.conn.commit()

		for fpath in glob.iglob(os.path.join(rootdir_inbox, '*', 'photometry', '*')):
			logger.info("="*72)
			logger.info(fpath)

			# Find the uploadlog corresponding to this file:
			db.cursor.execute("SELECT logid FROM flows.uploadlog WHERE uploadpath=%s;", [os.path.relpath(fpath, rootdir_inbox)])
			row = db.cursor.fetchone()
			if row is not None:
				uploadlogid = row['logid']
			else:
				uploadlogid = None

			# Only accept FITS file, or already compressed FITS files:
			if not fpath.endswith('.zip'):
				if uploadlogid:
					db.cursor.execute("UPDATE flows.uploadlog SET status='Invalid file type' WHERE logid=%s;", [uploadlogid])
					db.conn.commit()
				logger.error("Invalid file type: %s" % os.path.relpath(fpath, rootdir_inbox))
				continue

			# Get the name of the directory:
			# Not pretty, but it works...
			target_dirname = fpath[len(rootdir_inbox)+1:]
			target_dirname = target_dirname.split(os.path.sep)[0]

			# Convert directory name to target
			db.cursor.execute("SELECT targetid,target_name FROM flows.targets WHERE target_name=%s;", [target_dirname])
			row = db.cursor.fetchone()
			if row is None:
				logger.error('Could not find target: %s', target_dirname)
				continue
			targetid = row['targetid']
			targetname = row['target_name']


			try:
				with tempfile.TemporaryDirectory() as tmpdir:
					#
					tmpphotfile = os.path.join(tmpdir, 'photometry.ecsv')

					# Extract the ZIP file:
					with ZipFile(fpath, mode='r') as myzip:
						for member in myzip.infolist():
							# Remove any directory structure from the zip file:
							if member.is_dir():
								continue
							member.filename = os.path.basename(member.filename)

							# Due to security considerations, we only allow specific files
							# to be extracted:
							if member.filename == 'photometry.ecsv':
								myzip.extract(member, path=os.path.dirname(tmpdir))
							elif member.filename.endswith('.png') or member.filename.endswith('.log'):
								myzip.extract(member, path=os.path.dirname(tmpdir))

					# Check that the photometry ECSV file at least exists:
					if not os.path.isfile(tmpphotfile):
						raise FileNotFoundError("Photometry is not found")

					# Load photometry table:
					tab = Table.read(tmpphotfile, format='ascii.ecsv')
					fileid_img = int(tab.meta['fileid'])
					targetid_table = int(tab.meta['targetid'])

					assert targetid_table == targetid

					# Find out which version number to assign to file:
					db.cursor.execute("SELECT latest_version FROM flows.photometry_summary WHERE fileid_img=%s;", [fileid_img,])
					latest_version = db.cursor.fetchone()
					if latest_version is None:
						new_version = 1
					else:
						new_version = latest_version[0] + 1

					# Create a new path and filename that is slightly more descriptive:
					newpath = os.path.join(
						rootdir_archive,
						targetname,
						'{0:05d}'.format(fileid_img),
						'v{0:02d}'.format(new_version),
						'photometry-{0:s}-{1:05d}-v{2:02d}.ecsv'.format(targetname, fileid_img, new_version)
					)
					logger.info(newpath)

					if os.path.exists(newpath):
						print("Already exists")
						if uploadlogid:
							db.cursor.execute("UPDATE flows.uploadlog SET status='Already exists: file name' WHERE logid=%s;", [uploadlogid])
							db.conn.commit()
						continue

					archive, relpath = flows_get_archive_from_path(newpath, archives_list)

					db.cursor.execute("SELECT fileid FROM flows.files WHERE archive=%s AND path=%s;", [archive, relpath])
					if db.cursor.fetchone() is not None:
						print("ALREADY DONE")
						continue

					db.cursor.execute("SELECT * FROM flows.files WHERE fileid=%s;", [fileid_img])
					row = db.cursor.fetchone()
					site = row['site']

					assert targetid == row['targetid']
					assert tab.meta['photfilter'] == row['photfilter']

					# Optimize all the PNG files in the temp directory:
					for f in glob.iglob(os.path.join(tmpdir, '*.png')):
						optipng(f)

					# Copy the full directory to its new home:
					shutil.copytree(tmpdir, os.path.dirname(newpath))
					os.rename(os.path.join(os.path.dirname(newpath), 'photometry.ecsv'), newpath)

				# Get information about file:
				filesize = os.path.getsize(newpath)
				filehash = get_filehash(newpath)

				db.cursor.execute("INSERT INTO flows.files (archive,path,targetid,datatype,site,filesize,filehash,obstime,photfilter,version) VALUES (%(archive)s,%(relpath)s,%(targetid)s,%(datatype)s,%(site)s,%(filesize)s,%(filehash)s,%(obstime)s,%(photfilter)s,%(version)s) RETURNING fileid;", {
					'archive': archive,
					'relpath': relpath,
					'targetid': targetid,
					'datatype': 2,
					'site': site,
					'filesize': filesize,
					'filehash': filehash,
					'obstime': tab.meta['obstime-bmjd'],
					'photfilter': tab.meta['photfilter'],
					'version': new_version
				})
				fileid = db.cursor.fetchone()[0]

				# Add dependencies:
				db.cursor.execute("INSERT INTO flows.files_cross_assoc (fileid,associd) VALUES (%s,%s);", [fileid, fileid_img])
				if tab.meta['template'] is not None:
					db.cursor.execute("INSERT INTO flows.files_cross_assoc (fileid,associd) VALUES (%s,%s);", [fileid, tab.meta['template']])
				if tab.meta['diffimg'] is not None:
					db.cursor.execute("INSERT INTO flows.files_cross_assoc (fileid,associd) VALUES (%s,%s);", [fileid, tab.meta['diffimg']])

				indx_raw = (tab['starid'] == 0)
				indx_sub = (tab['starid'] == -1)
				print(tab[indx_raw])
				print(tab[indx_sub])

				phot_summary = {
					'fileid_img': fileid_img,
					'fileid_phot': fileid,
					'fileid_template': tab.meta['template'],
					'fileid_diffimg': tab.meta['diffimg'],
					'targetid': targetid,
					'obstime': tab.meta['obstime-bmjd'],
					'photfilter': tab.meta['photfilter'],
					'mag_raw': tab[indx_raw]['mag'],
					'mag_raw_error': tab[indx_raw]['mag_error'],
					'mag_sub': tab[indx_sub]['mag'],
					'mag_sub_error': tab[indx_sub]['mag_error'],
					'pipeline_version': tab.meta['version'],
					'latest_version': new_version
				}

				db.cursor.execute("SELECT * FROM flows.photometry_summary WHERE fileid_img=%s;", [fileid_img])
				if db.cursor.fetchone() is None:
					db.cursor.execute("INSERT INTO flows.photometry_summary (fileid_phot,fileid_img,fileid_template,fileid_diffimg,targetid,obstime,photfilter,mag,mag_error,pipeline_version,latest_version) VALUES (%(fileid_phot)s,%(fileid_img)s,%(fileid_template)s,%(fileid_diffimg)s,%(targetid)s,%(obstime)s,%(photfilter)s,%(mag)s,%(mag_error)s,%(pipeline_version)s,%(latest_version)s);", phot_summary)
				else:
					db.cursor.execute("UPDATE flows.photometry_summary SET fileid_phot=%(fileid_phot)s,targetid=%(targetid)s,fileid_template=%(fileid_template)s,fileid_diffimg=%(fileid_diffimg)s,obstime=%(obstime)s,photfilter=%(photfilter)s,mag=%(mag)s,mag_error=%(mag_error)s,pipeline_version=%(pipeline_version)s,latest_version=%(latest_version)s WHERE fileid_img=%(fileid_img)s;", phot_summary)

				if uploadlogid:
					db.cursor.execute("UPDATE flows.uploadlog SET fileid=%s,status='ok' WHERE logid=%s;", [fileid, uploadlogid])

				db.conn.commit()

			except:
				db.conn.rollback()
				if os.path.isdir(os.path.dirname(newpath)):
					shutil.rmtree(os.path.dirname(newpath))
				raise
			else:
				logger.info("DELETE THE ORIGINAL FILE")
				if os.path.isfile(fpath):
					os.remove(fpath)
				if uploadlogid:
					db.cursor.execute("UPDATE flows.uploadlog SET uploadpath=NULL WHERE logid=%s;", [uploadlogid])
					db.conn.commit()

#--------------------------------------------------------------------------------------------------
def cleanup_inbox():
	"""
	Cleanup of inbox directory
	"""
	rootdir_inbox = '/flows/inbox'

	for inputtype in ('science', 'templates', 'subtracted', 'photometry'):
		for dpath in glob.iglob(os.path.join(rootdir_inbox, '*', inputtype)):
			if not os.listdir(dpath):
				os.rmdir(dpath)

	for dpath in glob.iglob(os.path.join(rootdir_inbox, '*')):
		if os.path.isdir(dpath) and not os.listdir(dpath):
			os.rmdir(dpath)

	# Delete left-over files in the tables, that have been removed from disk:
	with AADC_DB() as db:
		db.cursor.execute("SELECT logid,uploadpath FROM flows.uploadlog WHERE uploadpath IS NOT NULL;")
		for row in db.cursor.fetchall():
			if not os.path.isfile(os.path.join(rootdir_inbox, row['uploadpath'])):
				print("DELETE THIS FILE: " + row['uploadpath'])
				db.cursor.execute("UPDATE flows.uploadlog SET uploadpath=NULL,status='File deleted' WHERE logid=%s;", [row['logid']])
				db.conn.commit()

#--------------------------------------------------------------------------------------------------
if __name__ == '__main__':

	logging_level = logging.INFO

	# Setup logging:
	formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
	console = logging.StreamHandler(sys.stdout)
	console.setFormatter(formatter)
	logger = logging.getLogger(__name__)
	if not logger.hasHandlers():
		logger.addHandler(console)
	logger.setLevel(logging_level)

	ingest_from_inbox()
	#ingest_photometry_from_inbox()
	cleanup_inbox()
