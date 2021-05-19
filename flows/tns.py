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
import astropy.units as u
import requests
import json
from bs4 import BeautifulSoup
import datetime
from .config import load_config

url_tns_api = 'https://www.wis-tns.org/api/get'
url_tns_search = 'https://www.wis-tns.org/search'

#--------------------------------------------------------------------------------------------------
class TNSConfigError(Exception):
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
		reply['internal_names'] = [name.strip() for name in reply['internal_names'].split(',') if name.strip()]
		return reply
	return None

#--------------------------------------------------------------------------------------------------
def tns_get_names(months=None, date_begin=None, date_end=None, z_min=None, z_max=None, objtype=None):
	"""Get SN names from TNS html table"""

	logger = logging.getLogger(__name__)

	# Change formats of input to be ready for query:
	if isinstance(date_begin, datetime.datetime):
		date_begin = date_begin.date().isoformat()
	elif isinstance(date_begin, datetime.date):
		date_end = date_begin.isoformat()
	if isinstance(date_end, datetime.datetime):
		date_end = date_end.date().isoformat()
	elif isinstance(date_end, datetime.date):
		date_end = date_end.isoformat()
	if isinstance(objtype, (list, tuple)):
		objtype = ','.join([str(o) for o in objtype])

	# Parameters for query:
	params = {
		'discovered_period_value': months,
		'discovered_period_units': 'months',
		'unclassified_at': 0,
		'classified_sne': 1,
		'include_frb': 0,
		#'name': ,
		'name_like': 0,
		'isTNS_AT': 'all',
		'public': 'all',
		#'ra':
		#'decl':
		#'radius':
		'coords_unit': 'arcsec',
		'reporting_groupid[]': 'null',
		'groupid[]': 'null',
		'classifier_groupid[]': 'null',
		'objtype[]': objtype,
		'at_type[]': 'null',
		'date_start[date]': date_begin,
		'date_end[date]': date_end,
		#'discovery_mag_min':
		#'discovery_mag_max':
		#'internal_name':
		#'discoverer':
		#'classifier':
		#'spectra_count':
		'redshift_min': z_min,
		'redshift_max': z_max,
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
		'display[redshift]': 1,
		'display[hostname]': 1,
		'display[host_redshift]': 1,
		'display[source_group_name]': 1,
		'display[classifying_source_group_name]': 1,
		'display[discovering_instrument_name]': 0,
		'display[classifing_instrument_name]': 0,
		'display[programs_name]': 0,
		'display[internal_name]': 1,
		'display[isTNS_AT]': 0,
		'display[public]': 1,
		'display[end_pop_period]': 0,
		'display[spectra_count]': 1,
		'display[discoverymag]': 1,
		'display[discmagfilter]': 1,
		'display[discoverydate]': 1,
		'display[discoverer]': 1,
		'display[remarks]': 0,
		'display[sources]': 0,
		'display[bibcode]': 0,
		'display[ext_catalogs]': 0
	}

	con = requests.get(url_tns_search, params=params)
	con.raise_for_status()

	# Parse output using BS4, find html table
	soup = BeautifulSoup(con.text, 'html.parser')
	tab = soup.find('table', 'results-table')
	if tab is None:
		logger.error('No HTML table obtained from query!')
		return None

	names = tab.find_all('td', 'cell-name')
	names_list = [name.text.replace(' ', '') for name in names if name.text.startswith('SN')]

	return names_list
