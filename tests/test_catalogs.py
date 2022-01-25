#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test API calls.

.. codeauthor:: Rasmus Handberg <rasmush@phys.au.dk>
"""

import pytest
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.table import Table
import conftest # noqa: F401
from flows import catalogs

#--------------------------------------------------------------------------------------------------
def test_query_simbad():

	# Coordinates around test-object (2019yvr):
	coo_centre = SkyCoord(
		ra=256.727512,
		dec=30.271482,
		unit='deg',
		frame='icrs'
	)

	results, simbad = catalogs.query_simbad(coo_centre)

	assert isinstance(results, Table)
	assert isinstance(simbad, SkyCoord)
	assert len(results) > 0
	results.pprint_all(50)

#--------------------------------------------------------------------------------------------------
def test_query_skymapper():

	# Coordinates around test-object (2021aess):
	coo_centre = SkyCoord(
		ra=53.4505,
		dec=-19.495725,
		unit='deg',
		frame='icrs'
	)

	results, skymapper = catalogs.query_skymapper(coo_centre)

	assert isinstance(results, Table)
	assert isinstance(skymapper, SkyCoord)
	assert len(results) > 0
	results.pprint_all(50)

#--------------------------------------------------------------------------------------------------

@pytest.mark.parametrize('ra,dec', [
	[256.727512, 30.271482], # 2019yvr
	[58.59512, -19.18172], # 2009D
])
def test_download_catalog(SETUP_CONFIG, ra, dec):

	# Check if CasJobs have been configured, and skip the entire test if it isn't.
	# This has to be done like this, to avoid problems when config.ini doesn't exist.
	try:
		catalogs.configure_casjobs()
	except catalogs.CasjobsError:
		pytest.skip("CasJobs not configured")

	# Coordinates around test-object (2019yvr):
	coo_centre = SkyCoord(
		ra=ra,
		dec=dec,
		unit='deg',
		frame='icrs'
	)

	tab = catalogs.query_all(coo_centre)
	print(tab)

	assert isinstance(tab, Table), "Should return a Table"
	results = catalogs.convert_table_to_dict(tab)

	assert isinstance(results, list), "Should return a list"
	for obj in results:
		assert isinstance(obj, dict), "Each element should be a dict"

		# Check columns:
		assert 'starid' in obj and obj['starid'] > 0
		assert 'ra' in obj and np.isfinite(obj['ra'])
		assert 'decl' in obj and np.isfinite(obj['decl'])
		assert 'pm_ra' in obj
		assert 'pm_dec' in obj
		assert 'gaia_mag' in obj
		assert 'gaia_bp_mag' in obj
		assert 'gaia_rp_mag' in obj
		assert 'gaia_variability' in obj
		assert 'B_mag' in obj
		assert 'V_mag' in obj
		assert 'u_mag' in obj
		assert 'g_mag' in obj
		assert 'r_mag' in obj
		assert 'i_mag' in obj
		assert 'z_mag' in obj
		assert 'H_mag' in obj
		assert 'J_mag' in obj
		assert 'K_mag' in obj

		# All values should be finite number or None:
		for key, val in obj.items():
			if key not in ('starid', 'gaia_variability'):
				assert val is None or np.isfinite(val), f"{key} is not a valid value: {val}"

	# TODO: Manually check a target from this position if the merge is correct.
	#assert False

#--------------------------------------------------------------------------------------------------
if __name__ == '__main__':
	pytest.main([__file__])
