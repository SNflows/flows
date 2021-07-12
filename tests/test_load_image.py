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
import os.path
import conftest # noqa: F401
from flows.load_image import load_image

#--------------------------------------------------------------------------------------------------
@pytest.mark.parametrize('fpath', ['SN2020aatc_K_20201213_495s.fits.gz'])
def test_load_image(fpath):

	# The test input directory containing the test-images:
	INPUT_DIR = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'input')

	# Load the image from the test-set:
	img = load_image(os.path.join(INPUT_DIR, fpath))

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
	assert isinstance(img.photfilter, str)
	assert isinstance(img.wcs, WCS)

#--------------------------------------------------------------------------------------------------
if __name__ == '__main__':
	pytest.main([__file__])
