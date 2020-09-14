#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

.. codeauthor:: Rasmus Handberg <rasmush@phys.au.dk>
"""

import astropy.units as u
from astropy.table import Table
from astropy.time import Time
import requests
from functools import lru_cache
from ..config import load_config

#--------------------------------------------------------------------------------------------------
@lru_cache(maxsize=10)
def get_catalog(target, radius=None, output='table'):
	"""

	Parameters:
		target (int or str):
		radius (float, optional): Radius around target in degrees to return targets for.
		outout (str, optional): Desired output format. Choises are 'table', 'dict', 'json'.
			Default='table'.

	Returns:
		dict: Dictionary with three members:
			- 'target': Information about target.
			- 'references': Table with information about reference stars close to target.
			- 'avoid': Table with stars close to target which should be avoided in FOV selection.

	.. codeauthor:: Rasmus Handberg <rasmush@phys.au.dk>
	"""

	assert output in ('table', 'json', 'dict'), "Invalid output format"

	# Get API token from config file:
	config = load_config()
	token = config.get('api', 'token', fallback=None)
	if token is None:
		raise Exception("No API token has been defined")

	#
	r = requests.get('https://flows.phys.au.dk/api/reference_stars.php',
		params={'target': target},
		headers={'Authorization': 'Bearer ' + token})
	r.raise_for_status()
	jsn = r.json()

	# Convert timestamps to actual Time objects:
	jsn['target']['inserted'] = Time(jsn['target']['inserted'], scale='utc')
	jsn['target']['discovery_date'] = Time(jsn['target']['discovery_date'], scale='utc')

	if output in ('json', 'dict'):
		return jsn

	dict_tables = {}

	tab = Table(
		names=('targetid', 'target_name', 'ra', 'decl', 'redshift', 'redshift_error', 'discovery_mag', 'catalog_downloaded', 'pointing_model_created', 'inserted', 'discovery_date', 'project', 'host_galaxy', 'ztf_id'),
		dtype=('int32', 'str', 'float64', 'float64', 'float32', 'float32', 'float32', 'bool', 'bool', 'object', 'object', 'str', 'str', 'str'),
		rows=[jsn['target']])

	tab['ra'].description = 'Right ascension'
	tab['ra'].unit = u.deg
	tab['decl'].description = 'Declination'
	tab['decl'].unit = u.deg
	dict_tables['target'] = tab

	for table_name in ('references', 'avoid'):
		tab = Table(
			names=('starid', 'ra', 'decl', 'pm_ra', 'pm_dec', 'gaia_mag', 'gaia_bp_mag', 'gaia_rp_mag', 'gaia_variability', 'B_mag', 'V_mag', 'H_mag', 'J_mag', 'K_mag', 'u_mag', 'g_mag', 'r_mag', 'i_mag', 'z_mag', 'distance'),
			dtype=('int64', 'float64', 'float64', 'float32', 'float32', 'float32', 'float32', 'float32', 'int32', 'float32', 'float32', 'float32', 'float32', 'float32', 'float32', 'float32', 'float32', 'float32', 'float32', 'float64'),
			rows=jsn[table_name])

		tab['starid'].description = 'Unique identifier in REFCAT2 catalog'
		tab['ra'].description = 'Right ascension'
		tab['ra'].unit = u.deg
		tab['decl'].description = 'Declination'
		tab['decl'].unit = u.deg
		tab['pm_ra'].description = 'Proper motion in right ascension'
		tab['pm_ra'].unit = u.mas/u.yr
		tab['pm_dec'].description = 'Proper motion in declination'
		tab['pm_dec'].unit = u.mas/u.yr
		tab['distance'].description = 'Distance from object to target'
		tab['distance'].unit = u.deg

		tab['gaia_mag'].description = 'Gaia G magnitude'
		tab['gaia_bp_mag'].description = 'Gaia Bp magnitude'
		tab['gaia_rp_mag'].description = 'Gaia Rp magnitude'
		tab['gaia_variability'].description = 'Gaia variability classification'
		tab['B_mag'].description = 'Johnson B magnitude'
		tab['V_mag'].description = 'Johnson V magnitude'
		tab['H_mag'].description = '2MASS H magnitude'
		tab['J_mag'].description = '2MASS J magnitude'
		tab['K_mag'].description = '2MASS K magnitude'
		tab['u_mag'].description = 'u magnitude'
		tab['g_mag'].description = 'g magnitude'
		tab['r_mag'].description = 'r magnitude'
		tab['i_mag'].description = 'i magnitude'
		tab['z_mag'].description = 'z magnitude'

		# Add some meta-data to the table as well:
		tab.meta['targetid'] = int(dict_tables['target']['targetid'])

		dict_tables[table_name] = tab

	return dict_tables

#--------------------------------------------------------------------------------------------------
def get_catalog_missing():

	# Get API token from config file:
	config = load_config()
	token = config.get('api', 'token', fallback=None)
	if token is None:
		raise Exception("No API token has been defined")

	#
	r = requests.get('https://flows.phys.au.dk/api/catalog_missing.php',
		headers={'Authorization': 'Bearer ' + token})
	r.raise_for_status()
	return r.json()
