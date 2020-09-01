# -*- coding: utf-8 -*-
"""

.. codeauthor:: Rasmus Handberg <rasmush@phys.au.dk>
"""

import os.path
import glob
from astropy.io import fits
from astropy.table import Table
import hashlib
import shutil
import getpass
import gzip
import re
from flows.aadc_db import AADC_DB
from flows.plots import plt, plot_image
from flows.load_image import load_image
from flows.config import load_config

#--------------------------------------------------------------------------------------------------
def flows_get_archive_from_path(fname, archives_list):

	archive = None
	relpath = None

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

	os.system('optipng -preserve "%s"' % output_fpath)

#--------------------------------------------------------------------------------------------------
def ingest_photometry(photfile=None):

	rootdir = '/flows/photometry'

	# Use local copy of archive if configured to do so:
	config = load_config()
	input_dir = config.get('photometry', 'output', fallback=None)
	if input_dir is None and photfile is None:
		raise ValueError()

	print(input_dir)

	if photfile is None:
		files = glob.glob(os.path.join(input_dir, '*', '*', 'photometry.ecsv'))
		print(files)
	elif os.path.isdir(photfile):
		files = glob.iglob(os.path.join(photfile, '**', 'photometry.ecsv'))
	elif os.path.isfile(photfile):
		files = [photfile]
	else:
		raise FileNotFoundError()

	passwd = getpass.getpass('Password for PostgreSQL: ')

	with AADC_DB(username='rasmush', password=passwd) as db:
		# Get list of archives:
		db.cursor.execute("SELECT archive,path FROM aadc.files_archives;")
		archives_list = db.cursor.fetchall()

		for fpath in files:
			tab = Table.read(fpath, format='ascii.ecsv')
			fileid_img = int(tab.meta['fileid'])
			targetid = int(tab.meta['targetid'])

			# Get target name from database:
			db.cursor.execute("SELECT target_name FROM flows.targets WHERE targetid=%s;", [targetid,])
			targetname = db.cursor.fetchone()[0]

			# Find out which version number to assign to file:
			db.cursor.execute("SELECT latest_version FROM flows.photometry_summary WHERE fileid_img=%s;", [fileid_img,])
			latest_version = db.cursor.fetchone()
			if latest_version is None:
				new_version = 1
			else:
				new_version = latest_version[0] + 1

			# Create a new path and filename that is slightly more descriptive:
			newpath = os.path.join(
				rootdir,
				targetname,
				'{0:05d}'.format(fileid_img),
				'v{0:02d}'.format(new_version),
				'photometry-{0:s}-{1:05d}-v{2:02d}.ecsv'.format(targetname, fileid_img, new_version)
			)
			archive, relpath = flows_get_archive_from_path(newpath, archives_list)

			print("="*72)
			print(fpath)
			print(newpath)

			db.cursor.execute("SELECT targetid FROM flows.files WHERE archive=%s AND path=%s;", [archive, relpath])
			if db.cursor.fetchone() is not None:
				print("ALREADY DONE")
				continue

			try:
				# Copy the full directory:
				shutil.copytree(os.path.dirname(fpath), os.path.dirname(newpath))
				os.remove(os.path.join(os.path.dirname(newpath), 'photometry.ecsv'))
				shutil.copy(fpath, newpath)

				filesize = os.path.getsize(newpath)
				filehash = get_filehash(newpath)

				db.cursor.execute("SELECT * FROM flows.files WHERE fileid=%s;", [fileid_img])
				row = db.cursor.fetchone()
				site = row['site']

				assert targetid == row['targetid']
				assert tab.meta['photfilter'] == row['photfilter']

				db.cursor.execute("INSERT INTO flows.files (archive,path,targetid,datatype,site,filesize,filehash,obstime,photfilter) VALUES (%(archive)s,%(relpath)s,%(targetid)s,%(datatype)s,%(site)s,%(filesize)s,%(filehash)s,%(obstime)s,%(photfilter)s) RETURNING fileid;", {
					'archive': archive,
					'relpath': relpath,
					'targetid': targetid,
					'datatype': 2,
					'site': site,
					'filesize': filesize,
					'filehash': filehash,
					'obstime': tab.meta['obstime-bmjd'],
					'photfilter': tab.meta['photfilter']
				})
				fileid = db.cursor.fetchone()[0]

				# Add dependencies:
				db.cursor.execute("INSERT INTO flows.files_cross_assoc (fileid,associd) VALUES (%s,%s);", [fileid, fileid_img])
				if tab.meta['template'] is not None:
					db.cursor.execute("INSERT INTO flows.files_cross_assoc (fileid,associd) VALUES (%s,%s);", [fileid, tab.meta['template']])

				phot_summary = {
					'fileid_img': fileid_img,
					'fileid_phot': fileid,
					'fileid_template': tab.meta['template'],
					'targetid': targetid,
					'obstime': tab.meta['obstime-bmjd'],
					'photfilter': tab.meta['photfilter'],
					'mag': tab[0]['mag'],
					'mag_error': tab[0]['mag_error'],
					'pipeline_version': tab.meta['version'],
					'latest_version': new_version
				}

				db.cursor.execute("SELECT * FROM flows.photometry_summary WHERE fileid_img=%s;", [fileid_img])
				if db.cursor.fetchone() is None:
					db.cursor.execute("INSERT INTO flows.photometry_summary (fileid_phot,fileid_img,fileid_template,targetid,obstime,photfilter,mag,mag_error,pipeline_version,latest_version) VALUES (%(fileid_phot)s,%(fileid_img)s,%(fileid_template)s,%(targetid)s,%(obstime)s,%(photfilter)s,%(mag)s,%(mag_error)s,%(pipeline_version)s,%(latest_version)s);", phot_summary)
				else:
					db.cursor.execute("UPDATE flows.photometry_summary SET fileid_phot=%(fileid_phot)s,targetid=%(targetid)s,fileid_template=%(fileid_template)s,obstime=%(obstime)s,photfilter=%(photfilter)s,mag=%(mag)s,mag_error=%(mag_error)s,pipeline_version=%(pipeline_version)s,latest_version=%(latest_version)s WHERE fileid_img=%(fileid_img)s;", phot_summary)

				db.conn.commit()

			except:
				db.conn.rollback()
				if os.path.isdir(os.path.dirname(newpath)):
					shutil.rmtree(os.path.dirname(newpath))
				raise
			else:
				print("DELETE THE ORIGINAL DIRECTORY")
				if os.path.isdir(os.path.dirname(fpath)):
					shutil.rmtree(os.path.dirname(fpath))

#---------------------------------------------------------------------------------------------------
def ingest_from_inbox():

	rootdir_inbox = '/flows/inbox'
	rootdir = '/flows/archive'

	# Check that root directories are available:
	if not os.path.isdir(rootdir_inbox):
		raise FileNotFoundError("INBOX does not exists")
	if not os.path.isdir(rootdir):
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

		for inputtype in ('science', 'templates', 'subtracted'):
			for fpath in glob.iglob(os.path.join(rootdir_inbox, '*', inputtype, '*')):
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
					print("Invalid file type: %s" % os.path.relpath(fpath, rootdir_inbox))
					continue

				# Get the name of the directory:
				# Not pretty, but it works...
				target_dirname = fpath[len(rootdir_inbox)+1:]
				target_dirname = target_dirname.split(os.path.sep)[0]

				# Convert directory name to target
				db.cursor.execute("SELECT targetid,target_name FROM flows.targets WHERE target_name=%s;", [target_dirname])
				row = db.cursor.fetchone()
				if row is None:
					raise Exception('Could not find target: %s', target_dirname)
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
						print("ORIGINAL SCIENCE IMAGE COULD NOT BE FOUND: %s" % os.path.basename(fpath))
						continue
					else:
						subtracted_original_fileid = subtracted_original_fileid[0]
				else:
					raise Exception("Not understood, Captain")

				print("="*72)
				print(fpath)
				print(newpath)

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
					print("Could not load FITS image")
					continue

				if img.site['siteid'] is None:
					raise Exception("Unknown SITE")

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
					print("DELETE THE ORIGINAL FILE")
					if os.path.isfile(newpath):
						os.remove(fpath)
					if uploadlogid:
						db.cursor.execute("UPDATE flows.uploadlog SET uploadpath=NULL WHERE logid=%s;", [uploadlogid])
						db.conn.commit()

	#----------------------------------------------------------------------------------------------
	# Cleanup of inbox directory:
	for inputtype in ('science', 'templates', 'subtracted'):
		for dpath in glob.iglob(os.path.join(rootdir_inbox, '*', inputtype)):
			if not os.listdir(dpath):
				os.rmdir(dpath)

	for dpath in glob.iglob(os.path.join(rootdir_inbox, '*')):
		if os.path.isdir(dpath) and not os.listdir(dpath):
			os.rmdir(dpath)

#--------------------------------------------------------------------------------------------------
if __name__ == '__main__':

	#ingest_photometry()
	ingest_from_inbox()
