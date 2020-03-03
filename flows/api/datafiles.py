# -*- coding: utf-8 -*-
"""

.. codeauthor:: Rasmus Handberg <rasmush@phys.au.dk>
"""

import requests
from datetime import datetime
from ..config import load_config

#--------------------------------------------------------------------------------------------------
def get_datafile(fileid):

	# Get API token from config file:
	config = load_config()
	token = config.get('api', 'token', fallback=None)
	if token is None:
		raise Exception("No API token has been defined")

	r = requests.get('https://flows.phys.au.dk/pipeline/datafiles.php',
		params={'fileid': fileid},
		headers={'Authorization': 'Bearer ' + token})
	r.raise_for_status()
	jsn = r.json()

	# Parse some of the fields to Python objects:
	jsn['inserted'] = datetime.strptime(jsn['inserted'], '%Y-%m-%d %H:%M:%S.%f')
	jsn['lastmodified'] = datetime.strptime(jsn['lastmodified'], '%Y-%m-%d %H:%M:%S.%f')

	return jsn

#--------------------------------------------------------------------------------------------------
def get_datafiles(targetid=None):

	# Get API token from config file:
	config = load_config()
	token = config.get('api', 'token', fallback=None)
	if token is None:
		raise Exception("No API token has been defined")

	r = requests.get('https://flows.phys.au.dk/pipeline/datafiles.php',
		params={'targetid': targetid},
		headers={'Authorization': 'Bearer ' + token})
	r.raise_for_status()
	jsn = r.json()

	return jsn
