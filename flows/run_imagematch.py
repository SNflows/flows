#!/usr/bin/env python
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
from setuptools import Distribution
from setuptools.command.install import install
from . import api

#--------------------------------------------------------------------------------------------------
class OnlyGetScriptPath(install):
    def run(self):
        # does not call install.run() by design
        self.distribution.install_scripts = self.install_scripts

def get_setuptools_script_dir():
    dist = Distribution({'cmdclass': {'install': OnlyGetScriptPath}})
    dist.dry_run = True  # not sure if necessary, but to be safe
    dist.parse_config_files()
    command = dist.get_command_obj('install')
    command.ensure_finalized()
    command.run()
    return dist.install_scripts

#--------------------------------------------------------------------------------------------------
def run_imagematch(datafile, target=None, star_coord=None, fwhm=None):
	"""

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
	imgmatch = os.path.join(get_setuptools_script_dir(), 'ImageMatch')
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
	if datafile['site'] == 1:
		config_file = os.path.join(__dir__, 'imagematch', 'imagematch_lcogt.cfg')
	elif datafile['site'] == 2:
		config_file = os.path.join(__dir__, 'imagematch', 'imagematch_hawki.cfg')
	else:
		config_file = os.path.join(__dir__, 'imagematch', 'imagematch_default.cfg')
	if not os.path.isfile(config_file):
		raise FileNotFoundError(config_file)

	# We will work in a temprary directory, since ImageMatch produces
	# a lot of extra output files that we don't want to have lying around
	# after it completes
	with tempfile.TemporaryDirectory() as tmpdir:

		# Copy the science and reference image to the temp dir:
		shutil.copy(reference_image, tmpdir)
		shutil.copy(science_image, tmpdir)

		# Construct the command to run ImageMatch:
		cmd = '"{python:s}" "{imgmatch:s}" -cfg "{config_file:s}" -snx {target_ra:.10f}d -sny {target_dec:.10f}d -m "{reference_image:s}" "{science_image:s}"'.format(
			python=sys.executable,
			imgmatch=imgmatch,
			config_file=config_file,
			reference_image=os.path.basename(reference_image),
			science_image=os.path.basename(science_image),
			target_ra=target['ra'],
			target_dec=target['decl']
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

		# Check the outputs from the subprocess:
		logger.info("Return code: %d", proc.returncode)
		logger.info("STDOUT:\n%s", stdout_data.strip())
		if stderr_data.strip() != '':
			logger.error("STDERR:\n%s", stderr_data.strip())
		if proc.returncode < 0:
			raise Exception("ImageMatch failed. Processed killed by OS with returncode %d." % proc.returncode)
		elif proc.returncode > 0:
			raise Exception("ImageMatch failed.")

		# Load the resulting difference image into memory:
		diffimg_path = os.path.join(tmpdir, os.path.splitext(os.path.basename(science_image))[0] + 'diff.fits')
		if not os.path.isfile(diffimg_path):
			raise FileNotFoundError(diffimg_path)

		with fits.open(diffimg_path, mode='readonly') as hdu:
			diffimg = np.asarray(hdu[0].data)

	return diffimg

