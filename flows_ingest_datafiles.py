# -*- coding: utf-8 -*-
"""

.. codeauthor:: Rasmus Handberg <rasmush@phys.au.dk>
"""

import os.path
import glob
from astropy.io import fits
from astropy.table import Table
from astropy.time import Time
import hashlib
import shutil
import getpass
from flows.aadc_db import AADC_DB
from flows.plots import plt, plot_image
from flows.load_image import load_image
from flows.config import load_config
from flows import api
import flows
import sys

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

	output_fpath = os.path.splitext(os.path.abspath(filepath))[0] + '.png'

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

	rootdir = '/aadc/flows/photometry'

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
		#db.cursor.execute("DELETE FROM flows.files WHERE datatype=2;")
		#db.conn.commit()

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

	targetid = 2
	targetname = 'SN2019yvr'

	rootdir_inbox = '/aadc/flows/inbox'
	rootdir = '/aadc/flows/archive'

	sites = api.sites.get_all_sites()
	site_keywords = {s['site_keyword']: s['siteid'] for s in sites}

	passwd = getpass.getpass('Password for PostgreSQL: ')

	with AADC_DB(username='rasmush', password=passwd) as db:

		# Get list of archives:
		db.cursor.execute("SELECT archive,path FROM aadc.files_archives;")
		archives_list = db.cursor.fetchall()

		for inputtype in ('observations', 'templates'):
			for fpath in glob.iglob(os.path.join(rootdir_inbox, '2019yvr', inputtype, '*.fits')):
				# Get the name of the directory:
				# Not pretty, but it works...
				target_dirname = fpath[len(rootdir_inbox):]
				target_dirname = target_dirname.split(os.path.sep)[0]
				print(target_dirname)

				# TODO: Convert directory name to target
				#db.cursor.execute("SELECT targetid FROM flows.targets WHERE target_name=%s;", [])
				#row = db.cursor.fetchone()
				#targetid = row['targetid']
				#targetname = row['target_name']

				if inputtype == 'observations':
					newpath = os.path.join(rootdir, targetname, os.path.basename(fpath))
					datatype = 1
				else:
					newpath = os.path.join(rootdir, targetname, inputtype, os.path.basename(fpath))
					datatype = 3

				print("="*72)
				print(fpath)
				print(newpath)

				if os.path.exists(newpath):
					print("Already exists")
					continue

				archive, relpath = flows_get_archive_from_path(newpath, archives_list)

				db.cursor.execute("SELECT fileid FROM flows.files WHERE archive=%s AND path=%s;", [archive, relpath])
				if db.cursor.fetchone() is not None:
					print("ALREADY DONE")
					continue

				filehash = get_filehash(fpath)

				db.cursor.execute("SELECT fileid FROM flows.files WHERE filehash=%s;", [filehash])
				if db.cursor.fetchone() is not None:
					print("ALREADY DONE: Filehash")
					continue

				try:
					# Copy the file to its new home:
					shutil.copy(fpath, newpath)

					filesize = os.path.getsize(fpath)

					with fits.open(newpath, mode='readonly', memmap=True) as hdu:
						hdr = hdu[0].header

						origin = hdr.get('ORIGIN')
						siteid = None
						if origin == 'LCOGT':
							site = hdr['SITE']
							siteid = site_keywords.get(site, None)

							#if os.path.basename(fpath).endswith('-e00.fits'):
							#	datatype = 0
							#else:
							#	datatype = 1

							#target = hdr['OBJECT']
							#db.cursor.execute("SELECT targetid FROM flows.targets WHERE target_name=%s;", [target])
							#row = db.cursor.fetchone()
							#if row is None:
							#	raise Exception()
							#targetid = row['targetid']

							photfilter = hdr['FILTER']
							obstime = hdr['MJD-OBS']
							exptime = hdr['EXPTIME']
						elif hdr.get('TELESCOP') == 'NOT' and hdr.get('INSTRUME') in ('ALFOSC FASU', 'ALFOSC_FASU') and hdr.get('OBS_MODE') == 'IMAGING':
							siteid = 5
							photfilter = {
								'B Bes': 'B',
								'V Bes': 'V',
								'g SDSS': 'gp',
								'r SDSS': 'rp',
								'i SDSS': 'ip',
								'u SDSS': 'up'
							}.get(hdr['FILTER'].replace('_', ' '), hdr['FILTER'])
							obstime = Time(hdr['DATE-AVG'], format='isot', scale='utc').mjd
							exptime = hdr['EXPTIME']
						elif hdr.get('FPA.TELESCOPE') == 'PS1' and hdr.get('FPA.INSTRUMENT') == 'GPC1':
							siteid = 7
							photfilter = {
								'g.00000': 'gp',
								'r.00000': 'rp',
								'i.00000': 'ip'
							}.get(hdr['FPA.FILTER'], hdr['FPA.FILTER'])
							obstime = hdr['MJD-OBS']
							exptime = hdr['EXPTIME']

					if siteid is None:
						print(hdr)
						raise Exception("Unknown SITE")

					if not fpath.endswith('-e00.fits'):
						create_plot(newpath)

					db.cursor.execute("INSERT INTO flows.files (archive,path,targetid,datatype,site,filesize,filehash,obstime,photfilter,exptime,available) VALUES (%(archive)s,%(relpath)s,%(targetid)s,%(datatype)s,%(site)s,%(filesize)s,%(filehash)s,%(obstime)s,%(photfilter)s,%(exptime)s,1);", {
						'archive': archive,
						'relpath': relpath,
						'targetid': targetid,
						'datatype': datatype,
						'site': siteid,
						'filesize': filesize,
						'filehash': filehash,
						'obstime': obstime,
						'photfilter': photfilter,
						'exptime': exptime
					})
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

#--------------------------------------------------------------------------------------------------
if __name__ == '__main__':

	ingest_photometry()
	#ingest_from_inbox()
