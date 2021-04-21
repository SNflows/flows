# -*- coding: utf-8 -*-
"""

.. codeauthor:: Rasmus Handberg <rasmush@phys.au.dk>
"""

import requests
from functools import lru_cache
import astropy.units as u
from astropy.coordinates import EarthLocation
from ..config import load_config

#--------------------------------------------------------------------------------------------------
@lru_cache(maxsize=10)
def get_site(siteid):

	# Get API token from config file:
	config = load_config()
	token = config.get('api', 'token', fallback=None)
	if token is None:
		raise Exception("No API token has been defined")

	r = requests.get('https://flows.phys.au.dk/api/sites.php',
		params={'siteid': siteid},
		headers={'Authorization': 'Bearer ' + token})
	r.raise_for_status()
	jsn = r.json()

	# Special derived objects:
	jsn['EarthLocation'] = EarthLocation(lat=jsn['latitude']*u.deg, lon=jsn['longitude']*u.deg, height=jsn['elevation']*u.m)

	return jsn

#--------------------------------------------------------------------------------------------------
@lru_cache(maxsize=1)
def get_all_sites():

	# Get API token from config file:
	config = load_config()
	token = config.get('api', 'token', fallback=None)
	if token is None:
		raise Exception("No API token has been defined")

	r = requests.get('https://flows.phys.au.dk/api/sites.php',
		headers={'Authorization': 'Bearer ' + token})
	r.raise_for_status()
	jsn = r.json()

	# Special derived objects:
	for site in jsn:
		site['EarthLocation'] = EarthLocation(lat=site['latitude']*u.deg, lon=site['longitude']*u.deg, height=site['elevation']*u.m)

	return jsn
