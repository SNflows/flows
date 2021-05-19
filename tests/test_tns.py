#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test TNS queries which rely on the TNS API.

.. codeauthor:: Rasmus Handberg <rasmush@phys.au.dk>
"""

import pytest
from astropy.coordinates import SkyCoord
from conftest import capture_cli
from flows import tns

#--------------------------------------------------------------------------------------------------
def test_tns_search(SETUP_CONFIG):

	coo_centre = SkyCoord(
		ra=191.283890127,
		dec=-0.45909033652,
		unit='deg',
		frame='icrs'
	)
	res = tns.tns_search(coo_centre)

	print(res)
	assert res[0]['objname'] == '2019yvr'
	assert res[0]['prefix'] == 'SN'

#--------------------------------------------------------------------------------------------------
def test_tns_get_obj(SETUP_CONFIG):

	res = tns.tns_get_obj('2019yvr')

	print(res)
	assert res['objname'] == '2019yvr'
	assert res['name_prefix'] == 'SN'

#--------------------------------------------------------------------------------------------------
def test_run_querytns(SETUP_CONFIG):

	# Run the command line interface:
	out, err, exitcode = capture_cli('run_querytns.py')
	assert exitcode == 0

#--------------------------------------------------------------------------------------------------
if __name__ == '__main__':
	pytest.main([__file__])
