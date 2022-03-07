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
from astropy.time import Time
import tempfile
import os
from conftest import capture_cli
from tendrils.utils import ztf


# --------------------------------------------------------------------------------------------------
def test_ztf_id():
    coo_centre = SkyCoord(ra=191.283890127, dec=-0.45909033652, unit='deg', frame='icrs')
    ztfid = ztf.query_ztf_id(coo_centre)
    assert ztfid == 'ZTF20aabqkxs'

    # With the correct discovery date we should get the same result:
    dd = Time('2019-12-27 12:30:14', format='iso', scale='utc')
    ztfid = ztf.query_ztf_id(coo_centre, discovery_date=dd)
    assert ztfid == 'ZTF20aabqkxs'

    # With a wrong discovery date, we should not get a ZTF id:
    dd = Time('2021-12-24 18:00:00', format='iso', scale='utc')
    ztfid = ztf.query_ztf_id(coo_centre, discovery_date=dd)
    assert ztfid is None

    coo_centre = SkyCoord(ra=181.6874198, dec=67.1649528, unit='deg', frame='icrs')
    ztfid = ztf.query_ztf_id(coo_centre)
    assert ztfid == 'ZTF21aatyplr'


# --------------------------------------------------------------------------------------------------
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


# --------------------------------------------------------------------------------------------------
@pytest.mark.parametrize('targetid', [2, 865])
def test_run_download_ztf(targetid):
    with tempfile.TemporaryDirectory() as tmpdir:
        # Nothing exists before running:
        assert len(os.listdir(tmpdir)) == 0

        # Run the command line interface:
        out, err, exitcode = capture_cli('run_download_ztf.py', [f'--target={targetid:d}', '-o', tmpdir])
        assert exitcode == 0

        # The output directory should now have two files:
        print(os.listdir(tmpdir))
        assert len(os.listdir(tmpdir)) == 2


# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    pytest.main([__file__])
