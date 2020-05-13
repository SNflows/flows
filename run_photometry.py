# -*- coding: utf-8 -*-
"""

.. codeauthor:: Rasmus Handberg <rasmush@phys.au.dk>
"""

import argparse
import logging
from flows import api, photometry, load_config

#--------------------------------------------------------------------------------------------------
if __name__ == '__main__':
	# Parse command line arguments:
	parser = argparse.ArgumentParser(description='Run photometry pipeline.')
	parser.add_argument('-d', '--debug', help='Print debug messages.', action='store_true')
	parser.add_argument('-q', '--quiet', help='Only report warnings and errors.', action='store_true')
	parser.add_argument('-o', '--overwrite', help='Overwrite existing results.', action='store_true')
	parser.add_argument('--fileid', type=int, default=None)
	parser.add_argument('--targetid', type=int, default=2)
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
	logger = logging.getLogger('flows')
	if not logger.hasHandlers():
		logger.addHandler(console)
	logger.setLevel(logging_level)

	if args.fileid is not None:
		# Run the specified fileid:
		fileids = [args.fileid]
	else:
		# Ask the API for a list of fileids which are yet to be processed:
		fileids = api.get_datafiles(targetid=args.targetid)
		print(fileids)

	config = load_config()
	output_folder_root = config.get('photometry', 'output', fallback='.')

	for fid in fileids:
		print("="*72)
		print(fid)
		print("="*72)

		datafile = api.get_datafile(fid)
		targetid = datafile['targetid']
		catalog = api.get_catalog(targetid, output='table')
		target = catalog['target'][0]
		target_name = str(target['target_name'])

		# Folder to save output:
		output_folder = os.path.join(output_folder_root, target_name, '%04d' % fid)
		os.makedirs(output_folder, exist_ok=True)

		# Also write any logging output to the
		_filehandler = logging.FileHandler(os.path.join(output_folder, 'photometry.log'), mode='w')
		_filehandler.setFormatter(formatter)
		_filehandler.setLevel(logging.INFO)
		logger.addHandler(_filehandler)

		try:
			photfile = photometry(fileid=fid, output_folder=output_folder)
		except:
			logger.exception("Photometry failed")
			photfile = None

		logger.removeHandler(_filehandler)

	# TODO: Ingest automatically?
