#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

.. codeauthor:: Rasmus Handberg <rasmush@phys.au.dk>
"""

import requests
from datetime import datetime
from functools import lru_cache
from ..config import load_config

#--------------------------------------------------------------------------------------------------
@lru_cache(maxsize=10)
def get_datafile(fileid):

	# Get API token from config file:
	config = load_config()
	token = config.get('api', 'token', fallback=None)
	if token is None:
		raise RuntimeError("No API token has been defined")

	r = requests.get('https://flows.phys.au.dk/api/datafiles.php',
		params={'fileid': fileid},
		headers={'Authorization': 'Bearer ' + token})
	r.raise_for_status()
	jsn = r.json()

	# Parse some of the fields to Python objects:
	jsn['inserted'] = datetime.strptime(jsn['inserted'], '%Y-%m-%d %H:%M:%S.%f')
	jsn['lastmodified'] = datetime.strptime(jsn['lastmodified'], '%Y-%m-%d %H:%M:%S.%f')

	return jsn

#--------------------------------------------------------------------------------------------------
def get_datafiles(targetid=None, filt=None):
	"""
	Get list of data file IDs to be processed.

	Parameters:
		targetid (int, optional): Target ID to process.
		filt (str, optional): Filter the returned list:
			- ``missing``: Return only data files that have not yet been processed.
			- ``'all'``: Return all data files.

	Returns:
		list: List of data files the can be processed.

	.. codeauthor:: Rasmus Handberg <rasmush@phys.au.dk>
	"""

	# Validate input:
	if filt is None:
		filt = 'missing'
	if filt not in ('missing', 'all', 'error'):
		raise ValueError("Invalid filter specified: '%s'" % filt)

	# Get API token from config file:
	config = load_config()
	token = config.get('api', 'token', fallback=None)
	if token is None:
		raise RuntimeError("No API token has been defined")

	params = {}
	if targetid is not None:
		params['targetid'] = targetid
	params['filter'] = filt

	r = requests.get('https://flows.phys.au.dk/api/datafiles.php',
		params=params,
		headers={'Authorization': 'Bearer ' + token})
	r.raise_for_status()
	jsn = r.json()

	return jsn
