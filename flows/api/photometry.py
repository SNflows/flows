#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Upload photometry results to Flows server.

.. codeauthor:: Rasmus Handberg <rasmush@phys.au.dk>
"""

import logging
import os
import zipfile
import requests
import tempfile
import shutil
import glob
from tqdm import tqdm
from astropy.table import Table
from .. import api
from ..config import load_config
from ..utilities import get_filehash

#--------------------------------------------------------------------------------------------------
def get_photometry(photid):
	"""
	Retrieve lightcurve from Flows server.

	Parameters:
		photid (int): Target to download lightcurve for.

	Returns:
		:class:`astropy.table.Table`: Table containing photometry.

	TODO:
		- Enable caching of files.

	.. codeauthor:: Rasmus Handberg <rasmush@phys.au.dk>
	"""

	# Get API token from config file:
	config = load_config()
	token = config.get('api', 'token', fallback=None)
	if token is None:
		raise RuntimeError("No API token has been defined")

	# Send query to the Flows API:
	params = {'fileid': photid}
	r = requests.get('https://flows.phys.au.dk/api/download_photometry.php',
		params=params,
		headers={'Authorization': 'Bearer ' + token})
	r.raise_for_status()

	# Create tempory directory and save the file into there,
	# then open the file as a Table:
	with tempfile.TemporaryDirectory() as tmpdir:
		tmpfile = os.path.join(tmpdir, f'photometry-{photid:d}.ecsv')
		with open(tmpfile, 'w') as fid:
			fid.write(r.text)

		tab = Table.read(tmpfile, format='ascii.ecsv')

	return tab


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
		raise RuntimeError("No API token has been defined")
	photdir_root = config.get('photometry', 'output', fallback='.')
	photdir = os.path.join(photdir_root, datafile['target_name'], '{0:05d}'.format(fileid))

	if not os.path.isdir(photdir):
		raise FileNotFoundError(photdir)

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

	# Path to the ZIP file:
	# TODO: Use tempfile instead?
	fname_zip = '{0:05d}.zip'.format(fileid)
	fpath_zip = os.path.join(photdir, fname_zip)

	try:
		# Create ZIP file with all the files:
		with zipfile.ZipFile(fpath_zip, 'w', allowZip64=True) as z:
			for f in tqdm(files, desc='Zipping {0:d}'.format(fileid), **tqdm_settings):
				logger.debug('Zipping %s', f)
				z.write(os.path.join(photdir, f), f)

		# Change the name of the uploaded file to contain the file hash:
		fname_zip = '{0:05d}-{1:s}.zip'.format(fileid, get_filehash(fpath_zip))

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
			raise RuntimeError("An error occurred while uploading photometry: " + r.text)
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
