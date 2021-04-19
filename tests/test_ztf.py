#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test ZTF queries which rely on the ZTF API.

.. codeauthor:: Rasmus Handberg <rasmush@phys.au.dk>
"""

import pytest
import numpy as np
from astropy.table import Table
from astropy.coordinates import SkyCoord
import conftest # noqa: F401
from flows import ztf

#--------------------------------------------------------------------------------------------------
def test_ztf_id():

	coo_centre = SkyCoord(
		ra=191.283890127,
		dec=-0.45909033652,
		unit='deg',
		frame='icrs'
	)
	ztfid = ztf.query_ztf_id(coo_centre)
	assert ztfid == 'ZTF20aabqkxs'

	coo_centre = SkyCoord(
		ra=181.6874198,
		dec=67.1649528,
		unit='deg',
		frame='icrs'
	)
	ztfid = ztf.query_ztf_id(coo_centre)
	assert ztfid == 'ZTF21aatyplr'

#--------------------------------------------------------------------------------------------------
@pytest.mark.parametrize('targetid', [2, 865])
def test_ztf_photometry(SETUP_CONFIG, targetid):

	tab = ztf.download_ztf_photometry(targetid)
	print(tab)

	assert isinstance(tab, Table)
	assert 'time' in tab.colnames
	assert 'photfilter' in tab.colnames
	assert 'mag' in tab.colnames
	assert 'mag_err' in tab.colnames
	assert np.all(np.isfinite(tab['time']))
	assert np.all(np.isfinite(tab['mag']))
	assert np.all(np.isfinite(tab['mag_err']))

#--------------------------------------------------------------------------------------------------
if __name__ == '__main__':
	pytest.main([__file__])
