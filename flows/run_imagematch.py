#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Flows photometry code.

.. codeauthor:: Rasmus Handberg <rasmush@phys.au.dk>
"""

import numpy as np
import logging
import os
import subprocess
import shlex
import tempfile
import sys
import shutil
import re
from astropy.io import fits
from astropy.wcs.utils import proj_plane_pixel_area
#from setuptools import Distribution
#from setuptools.command.install import install
from .load_image import load_image
from . import api

#--------------------------------------------------------------------------------------------------
#class OnlyGetScriptPath(install):
#	def run(self):
#		# does not call install.run() by design
#		self.distribution.install_scripts = self.install_scripts

#def get_setuptools_script_dir():
#	dist = Distribution({'cmdclass': {'install': OnlyGetScriptPath}})
#	dist.dry_run = True  # not sure if necessary, but to be safe
#	dist.parse_config_files()
#	command = dist.get_command_obj('install')
#	command.ensure_finalized()
#	command.run()
#	return dist.install_scripts

#--------------------------------------------------------------------------------------------------
def run_imagematch(datafile, target=None, star_coord=None, fwhm=None, pixel_scale=None):
	"""
	Run ImageMatch on a datafile.

	Parameters:
		datafile (dict): Data file to run ImageMatch on.
		target (:class:`astropy.table.Table`, optional): Target informaton.

	.. codeauthor:: Rasmus Handberg <rasmush@phys.au.dk>
	"""

	logger = logging.getLogger(__name__)

	if datafile.get('template') is None:
		raise ValueError("DATAFILE input does not specify a template to use.")

	# Extract paths to science and reference images:
	reference_image = os.path.join(datafile['archive_path'], datafile['template']['path'])
	science_image = os.path.join(datafile['archive_path'], datafile['path'])

	# If the target was not provided in the function call,
	# use the API to get the target information:
	if target is None:
		catalog = api.get_catalog(datafile['targetid'], output='table')
		target = catalog['target'][0]

	# Find the path to where the ImageMatch program is installed.
	# This is to avoid problems with it not being on the users PATH
	# and if the user is using some other version of the python executable.
	# TODO: There must be a better way of doing this!
	#imgmatch = os.path.join(get_setuptools_script_dir(), 'ImageMatch')
	if os.name == "nt":
		out = subprocess.check_output(["where", "ImageMatch"], universal_newlines=True)
		imgmatch = out.strip()
	else:
		out = subprocess.check_output(["whereis", "ImageMatch"], universal_newlines=True)
		out = re.match('ImageMatch: (.+)', out.strip())
		imgmatch = out.group(1)

	if not os.path.isfile(imgmatch):
		raise FileNotFoundError("ImageMatch not found")

	# Find the ImageMatch config file to use based on the site of the observations:
	__dir__ = os.path.dirname(os.path.abspath(__file__))
	if datafile['site'] in (1,3,4,6):
		config_file = os.path.join(__dir__, 'imagematch', 'imagematch_lcogt.cfg')
	elif datafile['site'] == 2:
		config_file = os.path.join(__dir__, 'imagematch', 'imagematch_hawki.cfg')
	elif datafile['site'] == 5:
		config_file = os.path.join(__dir__, 'imagematch', 'imagematch_alfosc.cfg')
	else:
		config_file = os.path.join(__dir__, 'imagematch', 'imagematch_default.cfg')
	if not os.path.isfile(config_file):
		raise FileNotFoundError(config_file)

	if pixel_scale is None:
		if datafile['site'] in (1,3,4,6):
			# LCOGT provides the pixel scale directly in the header
			pixel_scale = 'PIXSCALE'
		else:
			image = load_image(science_image)
			pixel_area = proj_plane_pixel_area(image.wcs)
			pixel_scale = np.sqrt(pixel_area)*3600 # arcsec/pixel
			logger.info("Calculated science image pixel scale: %f", pixel_scale)

	if datafile['template']['site'] in (1,3,4,6):
		# LCOGT provides the pixel scale directly in the header
		mscale = 'PIXSCALE'
	else:
		template = load_image(reference_image)
		template_pixel_area = proj_plane_pixel_area(template.wcs.celestial)
		mscale = np.sqrt(template_pixel_area)*3600 # arcsec/pixel
		logger.info("Calculated template pixel scale: %f", mscale)

	# Scale kernel radius with FWHM:
	if fwhm is None:
		kernel_radius = 9
	else:
		kernel_radius = max(9, int(np.ceil(1.5*fwhm)))
		if kernel_radius % 2 == 0:
			kernel_radius += 1

	# We will work in a temporary directory, since ImageMatch produces
	# a lot of extra output files that we don't want to have lying around
	# after it completes
	with tempfile.TemporaryDirectory() as tmpdir:

		# Copy the science and reference image to the temp dir:
		shutil.copy(reference_image, tmpdir)
		shutil.copy(science_image, tmpdir)

		# Construct the command to run ImageMatch:
		for match_threshold in (3.0, 5.0, 7.0, 10.0):
			cmd = '"{python:s}" "{imgmatch:s}" -cfg "{config_file:s}" -snx {target_ra:.10f}d -sny {target_dec:.10f}d -p {kernel_radius:d} -o {order:d} -s {match:f} -scale {pixel_scale:} -mscale {mscale:} -m "{reference_image:s}" "{science_image:s}"'.format(
				python=sys.executable,
				imgmatch=imgmatch,
				config_file=config_file,
				reference_image=os.path.basename(reference_image),
				science_image=os.path.basename(science_image),
				target_ra=target['ra'],
				target_dec=target['decl'],
				match=match_threshold,
				kernel_radius=kernel_radius,
				pixel_scale=pixel_scale,
				mscale=mscale,
				order=1
			)
			logger.info("Executing command: %s", cmd)

			# Run the command in a subprocess:
			cmd = shlex.split(cmd)
			proc = subprocess.Popen(cmd,
				cwd=tmpdir,
				stdout=subprocess.PIPE,
				stderr=subprocess.PIPE,
				universal_newlines=True)
			stdout_data, stderr_data = proc.communicate()
			returncode = proc.returncode
			proc.kill() # Cleanup - Is this really needed?

			# Check the outputs from the subprocess:
			logger.info("Return code: %d", returncode)
			logger.info("STDOUT:\n%s", stdout_data.strip())
			if stderr_data.strip() != '':
				logger.error("STDERR:\n%s", stderr_data.strip())
			if returncode < 0:
				raise Exception("ImageMatch failed. Processed killed by OS with returncode %d." % returncode)
			elif 'Failed object match... giving up.' in stdout_data:
				#raise Exception("ImageMatch giving up matching objects")
				continue
			elif returncode > 0:
				raise Exception("ImageMatch failed.")

			# Load the resulting difference image into memory:
			diffimg_name = re.sub(r'\.fits(\.gz|\.bz2)?$', r'diff.fits\1', os.path.basename(science_image))
			diffimg_path = os.path.join(tmpdir, diffimg_name)
			if not os.path.isfile(diffimg_path):
				raise FileNotFoundError(diffimg_path)

			break

		else:
			raise Exception("ImageMatch could not create difference image.")

		with fits.open(diffimg_path, mode='readonly') as hdu:
			diffimg = np.asarray(hdu[0].data)

	return diffimg
