# -*- coding: utf-8 -*-
"""
Created on Mon Feb 24 14:31:29 2020

@author: au195407
"""

import argparse
import logging
from flows import photometry

#--------------------------------------------------------------------------------------------------
if __name__ == '__main__':

	# Parse command line arguments:
	parser = argparse.ArgumentParser(description='Run TESS Photometry pipeline on single star.')
	parser.add_argument('-d', '--debug', help='Print debug messages.', action='store_true')
	parser.add_argument('-q', '--quiet', help='Only report warnings and errors.', action='store_true')
	parser.add_argument('-o', '--overwrite', help='Overwrite existing results.', action='store_true')
	parser.add_argument('-p', '--plot', help='Save plots when running.', action='store_true')
	args = parser.parse_args()

	# Set logging level:
	logging_level = logging.DEBUG
	if args.quiet:
		logging_level = logging.WARNING
	elif args.debug:
		logging_level = logging.DEBUG

	# Setup logging:
	formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
	console = logging.StreamHandler()
	console.setFormatter(formatter)
	logger = logging.getLogger('flows')
	logger.addHandler(console)
	logger.setLevel(logging_level)

	photometry.photometry(fileid=25)
