#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

.. codeauthor:: Rasmus Handberg <rasmush@phys.au.dk>
"""

import argparse
import logging
import os
import shutil
import functools
import multiprocessing
from flows import api, photometry, load_config

#--------------------------------------------------------------------------------------------------
def process_fileid(fid, output_folder_root=None, attempt_imagematch=True, autoupload=False):

	logger = logging.getLogger('flows')
	formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')

	datafile = api.get_datafile(fid)
	target_name = datafile['target_name']

	# Folder to save output:
	output_folder = os.path.join(output_folder_root, target_name, '%05d' % fid)

	photfile = None
	_filehandler = None
	try:
		# Set the status to indicate that we have started processing:
		api.set_photometry_status(fid, 'running')

		# Create the output directory if it doesn't exist:
		os.makedirs(output_folder, exist_ok=True)

		# Also write any logging output to the
		_filehandler = logging.FileHandler(os.path.join(output_folder, 'photometry.log'), mode='w')
		_filehandler.setFormatter(formatter)
		_filehandler.setLevel(logging.INFO)
		logger.addHandler(_filehandler)

		photfile = photometry(fileid=fid,
			output_folder=output_folder,
			attempt_imagematch=attempt_imagematch)

	except (SystemExit, KeyboardInterrupt):
		logger.error("Aborted by user or system.")
		if os.path.exists(output_folder):
			shutil.rmtree(output_folder, ignore_errors=True)
		photfile = None
		api.set_photometry_status(fid, 'abort')

	except: # noqa: E722, pragma: no cover
		logger.exception("Photometry failed")
		photfile = None
		api.set_photometry_status(fid, 'error')

	if _filehandler is not None:
		logger.removeHandler(_filehandler)

	if photfile is not None:
		if autoupload:
			api.upload_photometry(fid, delete_completed=True)
		api.set_photometry_status(fid, 'done')

	return photfile

#--------------------------------------------------------------------------------------------------
if __name__ == '__main__':
	# Parse command line arguments:
	parser = argparse.ArgumentParser(description='Run photometry pipeline.')
	parser.add_argument('-d', '--debug', help='Print debug messages.', action='store_true')
	parser.add_argument('-q', '--quiet', help='Only report warnings and errors.', action='store_true')
	parser.add_argument('-o', '--overwrite', help='Overwrite existing results.', action='store_true')

	group = parser.add_argument_group('Selecting which files to process')
	group.add_argument('--fileid', help="Process this file ID. Overrides all other filters.", type=int, default=None)
	group.add_argument('--targetid', help="Only process files from this target.", type=int, default=None)
	group.add_argument('--filter', type=str, default=None, choices=['missing','all'])

	group = parser.add_argument_group('Processing details')
	group.add_argument('--threads', type=int, default=1, help="Number of parallel threads to use.")
	group.add_argument('--no-imagematch', help="Disable ImageMatch.", action='store_true')
	group.add_argument('--autoupload', help="Automatically upload completed photometry to Flows website. Only do this, if you know what you are doing!", action='store_true')
	args = parser.parse_args()

	# Ensure that all input has been given:
	if not args.fileid and not args.targetid and args.filter is None:
		parser.error("Please select either a specific FILEID .")

	# Set logging level:
	logging_level = logging.INFO
	if args.quiet:
		logging_level = logging.WARNING
	elif args.debug:
		logging_level = logging.DEBUG

	# Number of threads to use:
	threads = args.threads
	if threads <= 0:
		threads = multiprocessing.cpu_count()

	# Setup logging:
	formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
	console = logging.StreamHandler()
	console.setFormatter(formatter)
	logger = logging.getLogger('flows')
	if not logger.hasHandlers():
		logger.addHandler(console)
	logger.setLevel(logging_level)

	if args.fileid is not None:
		# Run the specified fileid:
		fileids = [args.fileid]
	else:
		# Ask the API for a list of fileids which are yet to be processed:
		fileids = api.get_datafiles(targetid=args.targetid, filt=args.filter)

	config = load_config()
	output_folder_root = config.get('photometry', 'output', fallback='.')

	# Create function wrapper:
	process_fileid_wrapper = functools.partial(process_fileid,
		output_folder_root=output_folder_root,
		attempt_imagematch=not args.no_imagematch,
		autoupload=args.autoupload)

	if threads > 1:
		# Disable printing info messages from the parent function.
		# It is going to be all jumbled up anyway.
		#logger.setLevel(logging.WARNING)

		# There is more than one area to process, so let's start
		# a process pool and process them in parallel:
		with multiprocessing.Pool(threads) as pool:
			pool.map(process_fileid_wrapper, fileids)

	else:
		# Only single thread so simply run it directly:
		for fid in fileids:
			print("="*72)
			print(fid)
			print("="*72)
			process_fileid_wrapper(fid)
