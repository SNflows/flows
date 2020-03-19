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
#from . import api

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

	# Extract paths to science and reference images:
	reference_image = os.path.join(datafile['archive_path'], datafile['template'])
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
		tmpdir = os.path.abspath('./tmp12304543235')
		os.makedirs(tmpdir)

		# Copy the science and reference image to the temp dir:
		shutil.copy(reference_image, tmpdir)
		shutil.copy(science_image, tmpdir)

		# Construct the command to run ImageMatch:
		cmd = '"{python:s}" "{imgmatch:s}" -cfg "{config_file:s}" -v -snx {target_ra:.10f}d -sny {target_dec:.10f}d -m "{reference_image:s}" -sub "{science_image:s}"'.format(
			python=sys.executable,
			imgmatch=imgmatch,
			config_file=config_file,
			reference_image=os.path.basename(reference_image),
			science_image=os.path.basename(science_image),
			target_ra=target['ra'],
			target_dec=target['decl']
		)
		print(cmd)

		# Run the command in a subprocess:
		cmd = shlex.split(cmd)
		proc = subprocess.Popen(cmd,
			cwd=tmpdir,
			#stdout=subprocess.PIPE,
			stderr=subprocess.PIPE,
			universal_newlines=True)
		stdout_data, stderr_data = proc.communicate()

		# Check the outputs from the subprocess:
		print("Return code: %d" % proc.returncode)
		if proc.returncode < 0:
			raise Exception("")

		#if stderr_data is not None:
		print("ERROR: " + stderr_data)
		#raise Exception("ImageMatch failed")
		#print(stdout_data)

		print(os.listdir(tmpdir))

		# Load the resulting difference image into memory:
		diffimg_path = os.path.join(tmpdir, os.path.splitext(os.path.basename(science_image))[0] + 'diff.fits')
		if not os.path.isfile(diffimg_path):
			raise FileNotFoundError(diffimg_path)

		with fits.open(diffimg_path, mode='readonly') as hdu:
			diffimg = np.asarray(hdu[0].data)

	return diffimg

#--------------------------------------------------------------------------------------------------
if __name__ == '__main__':

	datafile = {'site': 1, 'archive_path': r'C:\Users\au195407\Documents\flows_archive', 'path': 'SN2019yvr/elp1m008-fa05-20200107-0439-e91.fits', 'template': 'SN2019yvr/templates/20170420T193430_Sloan_i_cpt1m010-fl16-20170420-0097-e91.fits'}
	target = {'ra': 191.283890127, 'decl': -0.45909033652}

	diffimg = run_imagematch(datafile, target)

	#ImageMatch -cfg "sample.cfg" *_g_*e91.fits -m "20170420T190451_Sloan_g_cpt1m010-fl16-20170420-0091-e91.fits"
	#ImageMatch -cfg "sample.cfg" *_r_*e91.fits -m "20170420T192210_Sloan_r_cpt1m010-fl16-20170420-0094-e91.fits"
	#ImageMatch -cfg "sample.cfg" *_i_*e91.fits -m "20170420T193430_Sloan_i_cpt1m010-fl16-20170420-0097-e91.fits"
