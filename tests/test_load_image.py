#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test loading of images.

.. codeauthor:: Rasmus Handberg <rasmush@phys.au.dk>
"""

import pytest
import numpy as np
from astropy.time import Time
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import os.path
import conftest # noqa: F401
from flows.api import get_filters
from flows.load_image import load_image

#--------------------------------------------------------------------------------------------------
@pytest.mark.parametrize('fpath,siteid', [
	['SN2020aatc_K_20201213_495s.fits.gz', 13],
	['ADP.2021-10-15T11_40_06.553.fits.gz', 2],
	#['TJO2459406.56826_V_imc.fits.gz', 22],
	#['lsc1m009-fa04-20210704-0044-e91_v1.fits.gz', 4],
	#['SN2021rcp_59409.931159242_B.fits.gz', 22],
	#['SN2021rhu_59465.86130221_B.fits.gz', 22],
	#['20200613_SN2020lao_u_stacked_meandiff.fits.gz', 1],
	#['2021aess_20220104_K.fits.gz', 5],
	#['2021aess_B01_20220207v1.fits.gz', 5],
])
def test_load_image(fpath, siteid):
	# Get list of all available filters:
	all_filters = set(get_filters().keys())

	# The test input directory containing the test-images:
	INPUT_DIR = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'input')

	# Target coordinates, only used for HAWKI image:
	target_coord = SkyCoord(
		ra=347.6230189,
		dec=7.5888196,
		unit='deg',
		frame='icrs')

	# Load the image from the test-set:
	img = load_image(os.path.join(INPUT_DIR, fpath), target_coord=target_coord)

	# Check the attributes of the image object:
	assert isinstance(img.image, np.ndarray)
	assert img.image.dtype in ('float32', 'float64')
	assert isinstance(img.mask, np.ndarray)
	assert img.mask.dtype == 'bool'
	assert isinstance(img.clean, np.ma.MaskedArray)
	assert img.clean.dtype == img.image.dtype
	assert isinstance(img.obstime, Time)
	assert isinstance(img.exptime, float)
	assert img.exptime > 0
	assert isinstance(img.wcs, WCS)
	assert isinstance(img.site, dict)
	assert img.site['siteid'] == siteid
	assert isinstance(img.photfilter, str)
	assert img.photfilter in all_filters

#--------------------------------------------------------------------------------------------------
if __name__ == '__main__':
	pytest.main([__file__])
