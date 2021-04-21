#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Query ZTF target information using ALeRCE API.
https://alerceapi.readthedocs.io/

.. codeauthor:: Emir Karamehmetoglu <emir.k@phys.au.dk>
.. codeauthor:: Rasmus Handberg <rasmush@phys.au.dk>
"""

import numpy as np
import astropy.units as u
from astropy.coordinates import Angle, SkyCoord
from astropy.table import Table
import datetime
import requests
from . import api

#--------------------------------------------------------------------------------------------------
def query_ztf_id(coo_centre, radius=1*u.arcmin):
	"""
	Query ALeRCE ZTF api to lookup ZTF identifier.

	In case multiple identifiers are found within the search cone, the one
	closest to the centre is returned.

	Parameters:
		coo_centre (:class:`astropy.coordinates.SkyCoord`): Coordinates of centre of search cone.
		radius (Angle, optional): Search radius. Default 1 arcmin.

	Returns:
		str: ZTF identifier.

	.. codeauthor:: Emir Karamehmetoglu <emir.k@phys.au.dk>
	.. codeauthor:: Rasmus Handberg <rasmush@phys.au.dk>
	"""

	if isinstance(radius, (float, int)):
		radius *= u.deg

	# Make json query for Alerce query API
	query = {
		'ra': coo_centre.ra.deg,
		'dec': coo_centre.dec.deg,
		'radius': Angle(radius).arcsec,
		'page_size': 20,
		'count': True
	}

	# Run http POST json query to alerce following their API
	res = requests.get('https://api.alerce.online/ztf/v1/objects', params=query)
	res.raise_for_status()
	jsn = res.json()

	# If nothing was found, return None:
	if jsn['total'] == 0:
		return None

	# Start by removing anything marked as likely stellar-like source:
	results = jsn['items']
	results = [itm for itm in results if not itm['stellar']]

	# Find target closest to the centre:
	coords = SkyCoord(
		ra=[itm['meanra'] for itm in results],
		dec=[itm['meandec'] for itm in results],
		unit='deg',
		frame='icrs')

	indx = np.argmin(coords.separation(coo_centre))

	return results[indx]['oid']

#--------------------------------------------------------------------------------------------------
def download_ztf_photometry(targetid):
	"""
	Download ZTF photometry from ALERCE API.

	Parameters:
		targetid (int): Target identifier.

	Returns:
		:class:`astropy.table.Table`: ZTF photometry table.

	.. codeauthor:: Emir Karamehmetoglu <emir.k@phys.au.dk>
	.. codeauthor:: Rasmus Handberg <rasmush@phys.au.dk>
	"""

	# Get target info from Flows API:
	tgt = api.get_target(targetid)
	oid = tgt['ztf_id']
	target_name = tgt['target_name']
	if oid is None:
		return None

	# Query ALERCE for detections of object based on oid
	res = requests.get(f'https://api.alerce.online/ztf/v1/objects/{oid:s}/detections')
	res.raise_for_status()
	jsn = res.json()

	# Create Astropy table, cut out the needed columns
	# and rename columns to something better for what we are doing:
	tab = Table(data=jsn)
	tab = tab[['fid', 'mjd', 'magpsf', 'sigmapsf']]
	tab.rename_column('fid', 'photfilter')
	tab.rename_column('mjd', 'time')
	tab.rename_column('magpsf', 'mag')
	tab.rename_column('sigmapsf', 'mag_err')

	# Remove bad values of time and magnitude:
	tab['time'] = np.asarray(tab['time'], dtype='float64')
	tab['mag'] = np.asarray(tab['mag'], dtype='float64')
	tab['mag_err'] = np.asarray(tab['mag_err'], dtype='float64')
	indx = np.isfinite(tab['time']) & np.isfinite(tab['mag']) & np.isfinite(tab['mag_err'])
	tab = tab[indx]

	# Replace photometric filter numbers with keywords used in Flows:
	photfilter_dict = {1: 'gp', 2: 'rp', 3: 'ip'}
	tab['photfilter'] = [photfilter_dict[fid] for fid in tab['photfilter']]

	# Sort the table on photfilter and time:
	tab.sort(['photfilter', 'time'])

	# Add meta information to table header:
	tab.meta['target_name'] = target_name
	tab.meta['targetid'] = targetid
	tab.meta['ztf_id'] = oid
	tab.meta['last_updated'] = datetime.datetime.now(tz=datetime.timezone.utc).isoformat()

	return tab
