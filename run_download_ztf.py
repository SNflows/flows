#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Download ZTF photometry from ALERCE API.

https://alerceapi.readthedocs.io/
"""

import argparse
import logging
import os.path
import numpy as np
from flows import ztf, api, load_config
from flows.plots import plt

#--------------------------------------------------------------------------------------------------
def main():
	# Parse command line arguments:
	parser = argparse.ArgumentParser(description='Download ZTF photometry.')
	parser.add_argument('-d', '--debug', help='Print debug messages.', action='store_true')
	parser.add_argument('-q', '--quiet', help='Only report warnings and errors.', action='store_true')
	parser.add_argument('-t', '--target', type=str, help='Target to download ZTF photometry for.', nargs='?', default=None)
	parser.add_argument('-o', '--output', type=str, help='Directory to save output to.', nargs='?', default=None)
	args = parser.parse_args()

	# Set logging level:
	logging_level = logging.INFO
	if args.quiet:
		logging_level = logging.WARNING
	elif args.debug:
		logging_level = logging.DEBUG

	# Setup logging:
	formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
	console = logging.StreamHandler()
	console.setFormatter(formatter)
	logger = logging.getLogger(__name__)
	if not logger.hasHandlers():
		logger.addHandler(console)
	logger.setLevel(logging_level)

	if args.output is None:
		config = load_config()
		output_dir = config.get('ztf', 'output_photometry', fallback='.')
	else:
		output_dir = args.output

	# Check that output directory exists:
	if not os.path.isdir(output_dir):
		parser.error("Output directory does not exist: '%s'" % output_dir)

	# Use API to get list of targets to process:
	if args.target is None:
		targets = api.get_targets()
	else:
		targets = [api.get_target(args.target)]

	# Loop through targets:
	for tgt in targets:
		logger.debug("Target: %s", tgt)

		# If there is no ZTF id, there is no need to try:
		if tgt['ztf_id'] is None:
			continue

		# Download ZTF photometry as Astropy Table:
		tab = ztf.download_ztf_photometry(tgt['targetid'])
		logger.debug("ZTF Photometry:\n%s", tab)
		if tab is None:
			continue

		# Write table to file:
		target_name = tab.meta['target_name']
		tab.write(os.path.join(output_dir, '{0:s}-ztf.ecsv'.format(target_name)),
			format='ascii.ecsv', delimiter=',')

		# Get LC data out and save as CSV files
		fig, ax = plt.subplots()
		for fid in np.unique(tab['photfilter']):
			band = tab[tab['photfilter'] == fid]
			ax.errorbar(band['mjd'], band['mag'], band['mag_err'],
				ls='none', marker='.', label=fid)

		ax.invert_yaxis()
		ax.set_title(target_name)
		ax.set_xlabel('Time (MJD)')
		ax.set_ylabel('Magnitude')
		ax.legend()
		fig.savefig(os.path.join(output_dir, '{0:s}-ztf.png'.format(target_name)), format='png')

#--------------------------------------------------------------------------------------------------
if __name__ == '__main__':
	main()
