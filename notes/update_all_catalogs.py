#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import logging
import sys
import os.path
from tqdm import tqdm
from astropy.coordinates import SkyCoord
if os.path.abspath('..') not in sys.path:
	sys.path.insert(0, os.path.abspath('..'))
import flows

#--------------------------------------------------------------------------------------------------
class TqdmLoggingHandler(logging.Handler):
	def __init__(self, *args, **kwargs):
		super().__init__(*args, **kwargs)

	def emit(self, record):
		try:
			msg = self.format(record)
			tqdm.tqdm.write(msg)
			self.flush()
		except (KeyboardInterrupt, SystemExit): # pragma: no cover
			raise
		except: # noqa: E722, pragma: no cover
			self.handleError(record)

#--------------------------------------------------------------------------------------------------
if __name__ == '__main__':

	# Setup logging:
	formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
	console = TqdmLoggingHandler()
	console.setFormatter(formatter)
	logger = logging.getLogger('flows')
	if not logger.hasHandlers():
		logger.addHandler(console)
	logger.setLevel(logging.INFO)

	for target in tqdm(flows.api.get_targets()):
		#if target['targetid'] != 1942:
		#	continue

		donefile = f"{target['targetid']:05d}.done"
		if not os.path.exists(donefile):
			flows.catalogs.download_catalog(target['targetid'], update_existing=True)

			open(donefile, 'w').close()
