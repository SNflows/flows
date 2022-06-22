"""
Test loading of images.
"""

import pytest
import numpy as np
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import os.path
import conftest  # noqa: F401
from tendrils import api
from flows.load_image import load_image

# Get list of all available filters:
ALL_FILTERS = set(api.get_filters().keys())


@pytest.mark.parametrize('fpath,siteid', [
    ['2020aatc/SN2020aatc_K_20201213_495s.fits.gz', 13],
    ['2021wyw/ADP.2021-10-15T11_40_06.553.fits.gz', 2], ])
def test_load_image(fpath, siteid):
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
    if img.mask is not None:
        assert isinstance(img.mask, np.ndarray)
        assert img.mask.dtype == 'bool'
    assert isinstance(img.clean, np.ma.MaskedArray)
    assert img.clean.dtype == img.image.dtype

    assert isinstance(img.exptime, float)
    assert img.exptime > 0.
    assert isinstance(img.wcs, WCS)

    assert isinstance(img.photfilter, str)
    assert img.photfilter in ALL_FILTERS


if __name__ == '__main__':
    pytest.main([__file__])
