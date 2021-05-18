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
import tempfile
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
	tqdm_settings = {'disable': None if logger.isEnabledFor(logging.INFO) else True}

	# Use API to get the datafile information:
	datafile = api.get_datafile(fileid)

	# Get API token from config file:
	config = load_config()
	token = config.get('api', 'token', fallback=None)
	if token is None:
		raise RuntimeError("No API token has been defined")
	photdir_root = config.get('photometry', 'output', fallback='.')

	# Find the photometry output directory for this fileid:
	photdir = os.path.join(photdir_root, datafile['target_name'], f'{fileid:05d}')
	if not os.path.isdir(photdir):
		# Do a last check, to ensure that we have not just added the wrong number of zeros
		# to the directory name:
		found_photdir = []
		for d in os.listdir(os.path.join(photdir_root, datafile['target_name'])):
			if d.isdigit() and int(d) == fileid and os.path.isdir(d):
				found_photdir.append(os.path.join(photdir_root, datafile['target_name'], d))
		# If we only found one, use it, otherwise throw an exception:
		if len(found_photdir) == 1:
			photdir = found_photdir[0]
		elif len(found_photdir) > 1:
			raise RuntimeError(f"Several photometry output found for fileid={fileid}. \
				You need to do a cleanup of the photometry output directories.")
		else:
			raise FileNotFoundError(photdir)

	# Make sure required files are actually there:
	photdir = os.path.abspath(photdir)
	files_existing = os.listdir(photdir)
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

	# Create ZIP file:
	with tempfile.TemporaryDirectory(prefix='flows-upload-') as tmpdir:
		# Create ZIP-file within the temp directory:
		fpath_zip = os.path.join(tmpdir, f'{fileid:05d}.zip')

		# Create ZIP file with all the files:
		with zipfile.ZipFile(fpath_zip, 'w', allowZip64=True) as z:
			for f in tqdm(files, desc=f'Zipping {fileid:d}', **tqdm_settings):
				logger.debug('Zipping %s', f)
				z.write(f, os.path.basename(f))

		# Change the name of the uploaded file to contain the file hash:
		fhash = get_filehash(fpath_zip)
		fname_zip = f'{fileid:05d}-{fhash:s}.zip'

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

	# If we have made it this far, the upload must have been a success:
	if delete_completed:
		if set([os.path.basename(f) for f in files]) == set(os.listdir(photdir)):
			logger.info("Deleting photometry from workdir: '%s'", photdir)
			shutil.rmtree(photdir, ignore_errors=True)
		else:
			logger.warning("Not deleting photometry from workdir: '%s'", photdir)
