#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

.. codeauthor:: Rasmus Handberg <rasmush@phys.au.dk>
"""

import requests
from functools import lru_cache
import astropy.units as u
from ..config import load_config

#--------------------------------------------------------------------------------------------------
@lru_cache(maxsize=10)
def get_filters():

	# Get API token from config file:
	config = load_config()
	token = config.get('api', 'token', fallback=None)
	if token is None:
		raise RuntimeError("No API token has been defined")

	r = requests.get('https://flows.phys.au.dk/api/filters.php',
		headers={'Authorization': 'Bearer ' + token})
	r.raise_for_status()
	jsn = r.json()

	# Add units:
	for f, val in jsn.items():
		if val.get('wavelength_center'):
			val['wavelength_center'] *= u.nm
		if val.get('wavelength_width'):
			val['wavelength_width'] *= u.nm

	return jsn
