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
from astropy.coordinates import Angle
from astropy.table import Table
import datetime
import requests
from . import api

#--------------------------------------------------------------------------------------------------
def query_ztf_id(coo_centre, radius=24*u.arcmin):
	"""
	Query ALeRCE ZTF api to lookup ZTF identifier.

	Parameters:
		coo_centre (:class:`astropy.coordinates.SkyCoord`): Coordinates of centre of search cone.
		radius (Angle, optional): Search radius. Default 24 arcmin.

	Returns:
		str: ZTF identifier.

	.. codeauthor:: Emir Karamehmetoglu <emir.k@phys.au.dk>
	.. codeauthor:: Rasmus Handberg <rasmush@phys.au.dk>
	"""

	if isinstance(radius, (float, int)):
		radius *= u.deg

	# Make json query for Alerce query API
	query = {
		"records_per_pages": 20,
		"query_parameters": {
			"coordinates": {
				"ra": coo_centre.ra.deg,
				"dec": coo_centre.dec.deg,
				"sr": Angle(radius).deg
			}
		}
	}

	# Run http POST json query to alerce following their API
	res = requests.post('https://ztf.alerce.online/query', json=query)
	res.raise_for_status()

	# If successful, get objectid of first object
	jsn = res.json()
	if jsn['total'] == 0:
		return None

	names = list(jsn['result'].keys())
	oid = names[0]
	return oid

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
	res = requests.post('https://ztf.alerce.online/get_detections', json={'oid': oid})
	res.raise_for_status()
	jsn = res.json()

	# Create Astropy table, cut out the needed columns
	# and rename columns to something better for what we are doing:
	tab = Table(data=jsn['result']['detections'])
	tab = tab[['fid', 'mjd', 'magpsf_corr', 'sigmapsf_corr']]
	tab.rename_column('fid', 'photfilter')
	tab.rename_column('mjd', 'time')
	tab.rename_column('magpsf_corr', 'mag')
	tab.rename_column('sigmapsf_corr', 'mag_err')

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
	tab.meta['last_updated'] = datetime.datetime.utcnow().isoformat()

	return tab
