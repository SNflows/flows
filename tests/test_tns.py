#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test TNS queries which rely on the TNS API.

.. codeauthor:: Rasmus Handberg <rasmush@phys.au.dk>
"""

import datetime
import logging
import os

import pytest
from astropy.coordinates import SkyCoord
from tendrils.utils import tns


# --------------------------------------------------------------------------------------------------
@pytest.mark.skipif(os.environ.get('CI') == 'true',
                    reason="Disabled on GitHub Actions to avoid too many requests HTTP error")
def test_tns_search():
    coo_centre = SkyCoord(ra=191.283890127, dec=-0.45909033652, unit='deg', frame='icrs')
    res = tns.tns_search(coo_centre)

    if res is None:
        raise ValueError("No results found.")
    assert res[0]['objname'] == '2019yvr'
    assert res[0]['prefix'] == 'SN'


# --------------------------------------------------------------------------------------------------
@pytest.mark.skipif(os.environ.get('CI') == 'true',
                    reason="Disabled on GitHub Actions to avoid too many requests HTTP error")
def test_tns_get_obj():
    res = tns.tns_get_obj('2019yvr')

    if res is None:
        raise ValueError("No results found.")

    assert res['objname'] == '2019yvr'
    assert res['name_prefix'] == 'SN'


# --------------------------------------------------------------------------------------------------
@pytest.mark.skipif(os.environ.get('CI') == 'true',
                    reason="Disabled on GitHub Actions to avoid too many requests HTTP error")
def test_tns_get_obj_noexist():
    res = tns.tns_get_obj('1892doesnotexist')
    assert res is None


# --------------------------------------------------------------------------------------------------
@pytest.mark.skipif(os.environ.get('CI') == 'true',
                    reason="Disabled on GitHub Actions to avoid too many requests HTTP error")
@pytest.mark.parametrize('date_begin,date_end',
                         [('2019-01-01', '2019-01-10'), (datetime.date(2019, 1, 1), datetime.date(2019, 1, 10)),
                          (datetime.datetime(2019, 1, 1, 12, 0), datetime.datetime(2019, 1, 10, 12, 0))])
def test_tns_getnames(date_begin, date_end):
    names = tns.tns_getnames(date_begin=date_begin, date_end=date_end, zmin=0, zmax=0.105, objtype=3)

    print(names)
    assert isinstance(names, list), "Should return a list"
    for n in names:
        assert isinstance(n, str), "Each element should be a string"
        assert n.startswith('SN'), "All names should begin with 'SN'"
    logging.debug(f"obtained names: {names}")
    assert 'SN2019A' in names, "SN2019A should be in the list"


# --------------------------------------------------------------------------------------------------
def test_tns_getnames_wronginput():
    # Wrong dates should result in ValueError:
    with pytest.raises(ValueError):
        tns.tns_getnames(date_begin=datetime.date(2019, 1, 1), date_end=datetime.date(2017, 1, 1))


# # --------------------------------------------------------------------------------------------------
# @pytest.mark.skipif(os.environ.get('CI') == 'true',
#                     reason="Disabled on GitHub Actions to avoid too many requests HTTP error")
# def test_run_querytns():
#     # Run the command line interface:
#     out, err, exitcode = capture_cli('run_querytns.py')
#     assert exitcode == 0
#
#
# # --------------------------------------------------------------------------------------------------
# if __name__ == '__main__':
#     pytest.main([__file__])
