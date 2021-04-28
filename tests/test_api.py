#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test API calls.

.. codeauthor:: Rasmus Handberg <rasmush@phys.au.dk>
"""

import pytest
import numpy as np
from astropy.coordinates import EarthLocation
from astropy.table import Table
import conftest # noqa: F401
from flows import api

#--------------------------------------------------------------------------------------------------
def test_api_get_targets(SETUP_CONFIG):

	tab = api.get_targets()
	print(tab)

	assert isinstance(tab, list)
	assert len(tab) > 0
	for target in tab:
		assert isinstance(target, dict)
		assert 'target_name' in target
		assert 'targetid' in target
		assert 'ra' in target
		assert 'decl' in target
		assert 'target_status' in target

#--------------------------------------------------------------------------------------------------
def test_api_get_target(SETUP_CONFIG):

	tab = api.get_target(2)
	print(tab)

	assert isinstance(tab, dict)
	assert tab['target_name'] == '2019yvr'
	assert tab['targetid'] == 2
	assert tab['target_status'] == 'target'
	assert tab['ztf_id'] == 'ZTF20aabqkxs'

#--------------------------------------------------------------------------------------------------
def test_api_get_datafiles(SETUP_CONFIG):

	tab = api.get_datafiles(targetid=2, filt='all')
	print(tab)
	assert isinstance(tab, list)
	assert len(tab) > 0
	for fid in tab:
		assert isinstance(fid, int)

	fileid = tab[0]
	tab = api.get_datafile(fileid)
	print(tab)
	assert tab['fileid'] == fileid
	assert tab['targetid'] == 2

#--------------------------------------------------------------------------------------------------
def test_api_get_filters(SETUP_CONFIG):

	tab = api.get_filters()
	print(tab)
	assert isinstance(tab, dict)
	for key, value in tab.items():
		assert isinstance(value, dict)
		assert value['photfilter'] == key
		assert 'wavelength_center' in value

#--------------------------------------------------------------------------------------------------
def test_api_get_sites(SETUP_CONFIG):

	tab = api.get_all_sites()
	print(tab)
	assert isinstance(tab, list)
	assert len(tab) > 0
	for site in tab:
		assert isinstance(site, dict)
		assert isinstance(site['siteid'], int)
		assert 'sitename' in site
		assert isinstance(site['EarthLocation'], EarthLocation)

	site0 = tab[0]
	print(site0)
	tab = api.get_site(site0['siteid'])
	print(tab)
	assert isinstance(tab, dict)
	assert tab == site0

#--------------------------------------------------------------------------------------------------
def test_api_get_catalog(SETUP_CONFIG):

	cat = api.get_catalog(2, output='table')
	print(cat)

	assert isinstance(cat, dict)

	target = cat['target']
	assert isinstance(target, Table)
	assert len(target) == 1
	assert target['targetid'] == 2
	assert target['target_name'] == '2019yvr'

	ref = cat['references']
	assert isinstance(ref, Table)

	avoid = cat['avoid']
	assert isinstance(avoid, Table)

#--------------------------------------------------------------------------------------------------
def test_api_get_lightcurve(SETUP_CONFIG):

	tab = api.get_lightcurve(2)
	print(tab)

	assert isinstance(tab, Table)
	assert len(tab) > 0
	assert 'time' in tab.colnames
	assert 'mag_raw' in tab.colnames

#--------------------------------------------------------------------------------------------------
def test_api_get_photometry(SETUP_CONFIG):

	tab = api.get_photometry(499)
	print(tab)

	# Basic tests of table:
	assert isinstance(tab, Table)
	assert len(tab) > 0
	assert 'starid' in tab.colnames
	assert 'ra' in tab.colnames
	assert 'decl' in tab.colnames
	assert 'mag' in tab.colnames
	assert 'mag_error' in tab.colnames
	assert np.sum(tab['starid'] == 0) == 1, "There should be one starid=0"

	# Meta-information:
	assert tab.meta['targetid'] == 2
	assert tab.meta['fileid'] == 179
	assert tab.meta['photfilter'] == 'B'

#--------------------------------------------------------------------------------------------------
if __name__ == '__main__':
	pytest.main([__file__])
