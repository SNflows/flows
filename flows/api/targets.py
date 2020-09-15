#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Get information about targets in Flows.

.. codeauthor:: Rasmus Handberg <rasmush@phys.au.dk>
"""

import requests
from functools import lru_cache
from ..config import load_config

#--------------------------------------------------------------------------------------------------
@lru_cache(maxsize=10)
def get_target(target):

	# Get API token from config file:
	config = load_config()
	token = config.get('api', 'token', fallback=None)
	if token is None:
		raise Exception("No API token has been defined")

	r = requests.get('https://flows.phys.au.dk/api/targets.php',
		params={'target': target},
		headers={'Authorization': 'Bearer ' + token})
	r.raise_for_status()
	jsn = r.json()

	return jsn

#--------------------------------------------------------------------------------------------------
@lru_cache(maxsize=1)
def get_targets():

	# Get API token from config file:
	config = load_config()
	token = config.get('api', 'token', fallback=None)
	if token is None:
		raise Exception("No API token has been defined")

	r = requests.get('https://flows.phys.au.dk/api/targets.php',
		headers={'Authorization': 'Bearer ' + token})
	r.raise_for_status()
	jsn = r.json()

	return jsn
