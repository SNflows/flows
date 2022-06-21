#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Tests of Photometry

.. codeauthor:: Rasmus Handberg <rasmush@phys.au.dk>
"""

import pytest
import conftest  # noqa: F401
from flows import photometry  # noqa: F401
import os
from pathlib import Path
from flows.fileio import DirectoriesDuringTest
from astropy.table import Table
from typing import Union


def test_import_photometry():
    pass


# fid = 9378  # NIR Cheap
# fid = 15005  # NIR Hawk-I Medium
# fid = 1446  # Optical and heavy.
def run_phot(fid):
    # The test input directory containing the test-images:
    input_dir = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'input/')
    output_dir = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'output/')
    tdir = DirectoriesDuringTest(input_dir, output_dir)
    photfile = photometry(fileid=fid, make_plots=False, directories=tdir)
    return photfile


def del_dir(target: Union[Path, str], only_if_empty: bool = False):
    """
    Delete a given directory and its subdirectories.

    :param target: The directory to delete
    :param only_if_empty: Raise RuntimeError if any file is found in the tree
    """
    target = Path(target).expanduser()
    assert target.is_dir()
    for p in sorted(target.glob('**/*'), reverse=True):
        if not p.exists():
            continue
        p.chmod(0o666)
        if p.is_dir():
            p.rmdir()
        else:
            if only_if_empty:
                raise RuntimeError(f'{p.parent} is not empty!')
            p.unlink()
    target.rmdir()


def clean_up_phot(filename):
    del_dir(Path(filename).parent)


@pytest.fixture(scope='session', autouse=True)
def photometry_fast():
    photfile = run_phot(fid=9378)
    yield photfile


@pytest.fixture(scope='session', autouse=True)
def phottable_fast(photometry_fast):
    yield Table.read(photometry_fast)
    clean_up_phot(photometry_fast)


def test_photometry_exists_fast(photometry_fast):
    assert os.path.isfile(photometry_fast)


def test_photometry_mag_fast(phottable_fast):
    sn_row = phottable_fast[0]
    assert sn_row['mag'] > 15.5 and sn_row['mag']< 15.9
    assert sn_row['mag_error'] > 0.11 and sn_row['mag_error'] < 0.18


# ------------------ SLOW

@pytest.mark.slow
@pytest.fixture(scope='session', autouse=True)
def photometry_slow_optical(fid=1446):
    photfile = run_phot(fid)
    yield photfile


@pytest.mark.slow
@pytest.fixture(scope='session', autouse=True)
def phottable_slow_optical(photometry_slow_optical):
    yield Table.read(photometry_slow_optical)
    clean_up_phot(photometry_slow_optical)


@pytest.mark.slow
def test_photometry_exists_slow(photometry_slow_optical):
    assert os.path.isfile(photometry_slow_optical)


@pytest.mark.slow
def test_photometry_mag_slow(phottable_slow_optical):
    sn_row = phottable_slow_optical[0]
    assert sn_row['mag'] > 17.35 and sn_row['mag'] < 17.5
    assert sn_row['mag_error'] > 0.003 and sn_row['mag_error'] < 0.02


if __name__ == '__main__':
    pytest.main([__file__])
