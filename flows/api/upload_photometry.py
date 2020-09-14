#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Upload photometry results to Flows server.

.. codeauthor:: Rasmus Handberg <rasmush@phys.au.dk>
"""

import logging
import os
import zipfile
import requests
import shutil
import glob
import time
from tqdm import tqdm
from .. import api
from ..config import load_config

#--------------------------------------------------------------------------------------------------
def upload_photometry(fileid, delete_completed=False):
	"""
	Upload photometry results to Flows server.

	This will make the uploaded photometry the active/newest/best photometry and
	be used in plots and shown on the website.

	Parameters:
		fileid (int): File ID of photometry to upload to server.
		delete_completed (bool, optional): Delete the photometry from the local
			working directory if the upload was successful. Default=False.

	.. codeauthor:: Rasmus Handberg <rasmush@phys.au.dk>
	"""

	logger = logging.getLogger(__name__)
	tqdm_settings = {'disable': not logger.isEnabledFor(logging.INFO)}

	# Use API to get the datafile information:
	datafile = api.get_datafile(fileid)

	# Get API token from config file:
	config = load_config()
	token = config.get('api', 'token', fallback=None)
	if token is None:
		raise Exception("No API token has been defined")
	photdir_root = config.get('photometry', 'output', fallback='.')
	photdir = os.path.join(photdir_root, datafile['target_name'], '{0:05d}'.format(fileid))

	if not os.path.isdir(photdir):
		raise FileNotFoundError(photdir)

	fname_zip = '{0:05d}.zip'.format(fileid)
	fpath_zip = os.path.join(photdir, fname_zip)
	files_existing = os.listdir(photdir)

	# Make sure required files are actually there:
	if 'photometry.ecsv' not in files_existing:
		raise FileNotFoundError(os.path.join(photdir, 'photometry.ecsv'))
	if 'photometry.log' not in files_existing:
		raise FileNotFoundError(os.path.join(photdir, 'photometry.log'))

	# Create list of files to be uploaded:
	files = [
		os.path.join(photdir, 'photometry.ecsv'),
		os.path.join(photdir, 'photometry.log')
	]
	files += glob.glob(os.path.join(photdir, '*.png'))

	try:
		# Create ZIP file with all the files:
		with zipfile.ZipFile(fpath_zip, 'w', allowZip64=True) as z:
			for f in tqdm(files, desc='Zipping {0:d}'.format(fileid), **tqdm_settings):
				logger.debug('Zipping %s', f)
				z.write(os.path.join(photdir, f), f)

		# Send file to the API:
		logger.info("Uploading to server...")
		with open(fpath_zip, 'rb') as fid:
			r = requests.post('https://flows.phys.au.dk/api/upload_photometry.php',
				params={'fileid': fileid},
				files={'file': (fname_zip, fid, 'application/zip')},
				headers={'Authorization': 'Bearer ' + token})

		# Check the returned data from the API:
		if r.text.strip() != 'OK':
			logger.error(r.text)
			raise Exception("An error occurred while uploading photometry: " + r.text)
		r.raise_for_status()

	finally:
		if os.path.isfile(fpath_zip):
			os.remove(fpath_zip)

	# If we have made it this far, the upload must have been a success:
	if delete_completed:
		if set([os.path.basename(f) for f in files]) == set(os.listdir(photdir)):
			logger.info("Deleting photometry from workdir: '%s'", photdir)
			shutil.rmtree(photdir, ignore_errors=True)
		else:
			logger.warning("Not deleting photometry from workdir: '%s'", photdir)
