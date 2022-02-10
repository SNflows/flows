#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
TNS API FUNCTIONS
Pre-provided helper functions for the TNS API
Obtained from https://wis-tns.weizmann.ac.il/content/tns-getting-started

.. codeauthor:: Emir Karamehmetoglu <emir.k@phys.au.dk>
.. codeauthor:: Rasmus Handberg <rasmush@phys.au.dk>
"""

import logging
from astropy.table import Table
import astropy.units as u
import requests
import json
import datetime
from .config import load_config

url_tns_api = 'https://www.wis-tns.org/api/get'
url_tns_search = 'https://www.wis-tns.org/search'

#--------------------------------------------------------------------------------------------------
class TNSConfigError(RuntimeError):
	pass

#--------------------------------------------------------------------------------------------------
def _load_tns_config():

	logger = logging.getLogger(__name__)

	config = load_config()
	api_key = config.get('TNS', 'api_key', fallback=None)
	if api_key is None:
		raise TNSConfigError("No TNS API-KEY has been defined in config")

	tns_bot_id = config.getint('TNS', 'bot_id', fallback=93222)
	tns_bot_name = config.get('TNS', 'bot_name', fallback='AUFLOWS_BOT')
	tns_user_id = config.getint('TNS', 'user_id', fallback=None)
	tns_user_name = config.get('TNS', 'user_name', fallback=None)

	if tns_user_id and tns_user_name:
		logger.debug('Using TNS credentials: user=%s', tns_user_name)
		user_agent = 'tns_marker{"tns_id":' + str(tns_user_id) + ',"type":"user","name":"' + tns_user_name + '"}'
	elif tns_bot_id and tns_bot_name:
		logger.debug('Using TNS credentials: bot=%s', tns_bot_name)
		user_agent = 'tns_marker{"tns_id":' + str(tns_bot_id) + ',"type":"bot","name":"' + tns_bot_name + '"}'
	else:
		raise TNSConfigError("No TNS bot_id or bot_name has been defined in config")

	return {
		'api-key': api_key,
		'user-agent': user_agent
	}

#--------------------------------------------------------------------------------------------------
def tns_search(coord=None, radius=3*u.arcsec, objname=None, internal_name=None):
	"""
	Cone-search TNS for object near coordinate.

	Parameters:
		coord (:class:`astropy.coordinates.SkyCoord`): Central coordinate to search around.
		radius (Angle, optional): Radius to search around ``coord``.
		objname (str, optional): Search on object name.
		internal_name (str, optional): Search on internal name.

	Returns:
		dict: Dictionary with TSN response.
	"""

	# API key for Bot
	tnsconf = _load_tns_config()

	# change json_list to json format
	json_file = {
		'radius': radius.to('arcsec').value,
		'units': 'arcsec',
		'objname': objname,
		'internal_name': internal_name
	}
	if coord:
		json_file['ra'] = coord.icrs.ra.deg
		json_file['dec'] = coord.icrs.dec.deg

	# construct the list of (key,value) pairs
	headers = {'user-agent': tnsconf['user-agent']}
	search_data = [
		('api_key', (None, tnsconf['api-key'])),
		('data', (None, json.dumps(json_file)))
	]

	# search obj using request module
	res = requests.post(url_tns_api + '/search', files=search_data, headers=headers)
	res.raise_for_status()
	parsed = res.json()
	data = parsed['data']

	if 'reply' in data:
		return data['reply']
	return None

#--------------------------------------------------------------------------------------------------
def tns_get_obj(name):
	"""
	Search TNS for object by name.

	Parameters:
		name (str): Object name to search for.

	Returns:
		dict: Dictionary with TSN response.
	"""

	# API key for Bot
	tnsconf = _load_tns_config()

	# construct the list of (key,value) pairs
	headers = {'user-agent': tnsconf['user-agent']}
	params = {'objname': name, 'photometry': '0', 'spectra': '0'}
	get_data = [
		('api_key', (None, tnsconf['api-key'])),
		('data', (None, json.dumps(params)))
	]

	# get obj using request module
	res = requests.post(url_tns_api + '/object', files=get_data, headers=headers)
	res.raise_for_status()
	parsed = res.json()
	data = parsed['data']

	if 'reply' in data:
		reply = data['reply']
		if not reply:
			return None
		if 'objname' not in reply: # Bit of a cheat, but it is simple and works
			return None

		reply['internal_names'] = [name.strip() for name in reply['internal_names'].split(',') if name.strip()]
		return reply
	return None

#--------------------------------------------------------------------------------------------------
def tns_getnames(months=None, date_begin=None, date_end=None, zmin=None, zmax=None, objtype=[3, 104]):
	"""
	Get SN names from TNS.

	Parameters:
		months (int, optional): Only return objects reported within the last X months.
		date_begin (date, optional): Discovery date begin.
		date_end (date, optional): Discovery date end.
		zmin (float, optional): Minimum redshift.
		zmax (float, optional): Maximum redshift.
		objtype (list, optional): Constraint object type.
			Default is to query for
			- 3: SN Ia
			- 104: SN Ia-91T-like

	Returns:
		list: List of names fulfilling search criteria.
	"""

	logger = logging.getLogger(__name__)

	# Change formats of input to be ready for query:
	if isinstance(date_begin, datetime.datetime):
		date_begin = date_begin.date()
	elif isinstance(date_begin, str):
		date_begin = datetime.datetime.strptime(date_begin, '%Y-%m-%d').date()

	if isinstance(date_end, datetime.datetime):
		date_end = date_end.date()
	elif isinstance(date_end, str):
		date_end = datetime.datetime.strptime(date_end, '%Y-%m-%d').date()

	if isinstance(objtype, (list, tuple)):
		objtype = ','.join([str(o) for o in objtype])

	# Do some sanity checks:
	if date_end < date_begin:
		raise ValueError("Dates are in the wrong order.")

	date_now = datetime.datetime.now(datetime.timezone.utc).date()
	if months is not None and date_end is not None and date_end < date_now - datetime.timedelta(days=months*30):
		logger.warning('Months limit restricts days_begin, consider increasing limit_months.')

	# API key for Bot
	tnsconf = _load_tns_config()

	# Parameters for query:
	params = {
		'discovered_period_value': months, # Reported Within The Last
		'discovered_period_units': 'months',
		'unclassified_at': 0, # Limit to unclasssified ATs
		'classified_sne': 1, # Limit to classified SNe
		'include_frb': 0, # Include FRBs
		#'name': ,
		'name_like': 0,
		'isTNS_AT': 'all',
		'public': 'all',
		#'ra':
		#'decl':
		#'radius':
		#'coords_unit': 'arcsec',
		'reporting_groupid[]': 'null',
		'groupid[]': 'null',
		'classifier_groupid[]': 'null',
		'objtype[]': objtype,
		'at_type[]': 'null',
		'date_start[date]': date_begin.isoformat(),
		'date_end[date]': date_end.isoformat(),
		#'discovery_mag_min':
		#'discovery_mag_max':
		#'internal_name':
		#'discoverer':
		#'classifier':
		#'spectra_count':
		'redshift_min': zmin,
		'redshift_max': zmax,
		#'hostname':
		#'ext_catid':
		#'ra_range_min':
		#'ra_range_max':
		#'decl_range_min':
		#'decl_range_max':
		'discovery_instrument[]': 'null',
		'classification_instrument[]': 'null',
		'associated_groups[]': 'null',
		#'at_rep_remarks':
		#'class_rep_remarks':
		#'frb_repeat': 'all'
		#'frb_repeater_of_objid':
		'frb_measured_redshift': 0,
		#'frb_dm_range_min':
		#'frb_dm_range_max':
		#'frb_rm_range_min':
		#'frb_rm_range_max':
		#'frb_snr_range_min':
		#'frb_snr_range_max':
		#'frb_flux_range_min':
		#'frb_flux_range_max':
		'num_page': 500,
		'display[redshift]': 0,
		'display[hostname]': 0,
		'display[host_redshift]': 0,
		'display[source_group_name]': 0,
		'display[classifying_source_group_name]': 0,
		'display[discovering_instrument_name]': 0,
		'display[classifing_instrument_name]': 0,
		'display[programs_name]': 0,
		'display[internal_name]': 0,
		'display[isTNS_AT]': 0,
		'display[public]': 0,
		'display[end_pop_period]': 0,
		'display[spectra_count]': 0,
		'display[discoverymag]': 0,
		'display[discmagfilter]': 0,
		'display[discoverydate]': 0,
		'display[discoverer]': 0,
		'display[remarks]': 0,
		'display[sources]': 0,
		'display[bibcode]': 0,
		'display[ext_catalogs]': 0,
		'format': 'csv'
	}

	# Query TNS for names:
	headers = {'user-agent': tnsconf['user-agent']}
	con = requests.get(url_tns_search, params=params, headers=headers)
	con.raise_for_status()

	# Parse the CSV table:
	# Ensure that there is a newline in table string.
	# AstroPy uses this to distinguish file-paths from pure-string inputs:
	text = str(con.text) + "\n"
	tab = Table.read(text,
		format='ascii.csv',
		guess=False,
		delimiter=',',
		quotechar='"',
		header_start=0,
		data_start=1)

	# Pull out the names only if they begin with "SN":
	names_list = [name.replace(' ', '') for name in tab['Name'] if name.startswith('SN')]
	names_list = sorted(names_list)

	return names_list
