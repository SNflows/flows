#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Get information about targets in Flows.

.. codeauthor:: Rasmus Handberg <rasmush@phys.au.dk>
"""

import re
from datetime import datetime
import pytz
from astropy.time import Time
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
		raise RuntimeError("No API token has been defined")

	r = requests.get('https://flows.phys.au.dk/api/targets.php',
		params={'target': target},
		headers={'Authorization': 'Bearer ' + token})
	r.raise_for_status()
	jsn = r.json()

	# Parse some of the fields to Python objects:
	jsn['inserted'] = datetime.strptime(jsn['inserted'], '%Y-%m-%d %H:%M:%S.%f')
	if jsn['discovery_date']:
		jsn['discovery_date'] = Time(jsn['discovery_date'], format='iso', scale='utc')

	return jsn

#--------------------------------------------------------------------------------------------------
@lru_cache(maxsize=1)
def get_targets():

	# Get API token from config file:
	config = load_config()
	token = config.get('api', 'token', fallback=None)
	if token is None:
		raise RuntimeError("No API token has been defined")

	r = requests.get('https://flows.phys.au.dk/api/targets.php',
		headers={'Authorization': 'Bearer ' + token})
	r.raise_for_status()
	jsn = r.json()

	# Parse some of the fields to Python objects:
	for tgt in jsn:
		tgt['inserted'] = datetime.strptime(tgt['inserted'], '%Y-%m-%d %H:%M:%S.%f')
		if tgt['discovery_date']:
			tgt['discovery_date'] = Time(tgt['discovery_date'], format='iso', scale='utc')

	return jsn

#--------------------------------------------------------------------------------------------------
def add_target(name, coord, redshift=None, redshift_error=None, discovery_date=None,
	discovery_mag=None, host_galaxy=None, ztf=None, status='candidate', project='flows'):
	"""
	Add new candidate or target.

	Coordinates are specified using an Astropy SkyCoord object, which can be
	created in the following way:

	coord = SkyCoord(ra=19.1, dec=89.00001, unit='deg', frame='icrs')

	The easiest way is to specify ``discovery_date`` as an Astropy Time object:

	discovery_date = Time('2020-01-02 00:00:00', format='iso', scale='utc')

	Alternatively, you can also specify it as a :class:`datetime.datetime` object,
	but some care has to be taken with specifying the correct timezone:

	discovery_date = datetime.strptime('2020-01-02 00:00:00', '%Y-%m-%d %H:%M:%S%f')
	discovery_date = pytz.timezone('America/New_York').localize(ddate)

	Lastly, it can be given as a simple date-string of the following form,
	but here the data has to be given in UTC:

	discovery_date = '2020-01-02 23:56:02.123'

	Parameters:
		name (str): Name of target. Must be of the form "YYYYxyz", where YYYY is the year,
			and xyz are letters.
		coord (:class:ʼastropy.coordinates.SkyCoordʼ): Sky coordinates of target.
		redshift (float, optional): Redshift.
		redshift_error (float, optional): Uncertainty on redshift.
		discovery_date (:class:`astropy.time.Time`, :class:`datetime.datetime` or str, optional):
		discovery_mag (float, optional): Magnitude at time of discovery.
		host_galaxy (str, optional): Host galaxy name.
		ztf (str, optional): ZTF identifier.
		status (str, optional):
		project (str, optional):

	Returns:
		int: New target identifier in Flows system.

	.. codeauthor:: Rasmus Handberg <rasmush@phys.au.dk>
	"""
	# Check and convert input:
	if not re.match(r'^[12]\d{3}([A-Z]|[a-z]{2,4})$', name.strip()):
		raise ValueError("Invalid target name.")

	if redshift is None and redshift_error is not None:
		raise ValueError("Redshift error specified without redshift value")

	if isinstance(discovery_date, Time):
		discovery_date = discovery_date.utc.iso
	elif isinstance(discovery_date, datetime):
		discovery_date = discovery_date.astimezone(pytz.timezone('UTC'))
		discovery_date = discovery_date.strftime('%Y-%m-%d %H:%M:%S%f')
	elif isinstance(discovery_date, str):
		discovery_date = datetime.strptime(discovery_date, '%Y-%m-%d %H:%M:%S%f')
		discovery_date = discovery_date.strftime('%Y-%m-%d %H:%M:%S%f')

	if status not in ('candidate', 'target'):
		raise ValueError("Invalid target status.")

	# Get API token from config file:
	config = load_config()
	token = config.get('api', 'token', fallback=None)
	if token is None:
		raise RuntimeError("No API token has been defined")

	# Gather parameters to be sent to API:
	params = {
		'targetid': 0,
		'target_name': name.strip(),
		'ra': coord.icrs.ra.deg,
		'decl': coord.icrs.dec.deg,
		'redshift': redshift,
		'redshift_error': redshift_error,
		'discovery_date': discovery_date,
		'discovery_mag': discovery_mag,
		'host_galaxy': host_galaxy,
		'project': project,
		'ztf_id': ztf,
		'target_status': status
	}

	# Post the request to the API:
	r = requests.post('https://flows.phys.au.dk/api/targets_add.php',
		data=params,
		headers={'Authorization': 'Bearer ' + token})
	r.raise_for_status()
	jsn = r.json()

	# Check for errors:
	if jsn['errors'] is not None:
		raise RuntimeError(f"Adding target '{name}' resulted in an error: {jsn['errors']}")

	return int(jsn['targetid'])
