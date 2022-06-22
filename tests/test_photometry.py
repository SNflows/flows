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
    results_table = photometry(fileid=fid, make_plots=False, directories=tdir)
    tdir.set_output_dirs(results_table.meta['target_name'], fid)
    return results_table, tdir


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
    result_table, tdir = run_phot(fid=9378)
    yield result_table, tdir.photometry_path
    clean_up_phot(tdir.photometry_path)


def test_photometry_exists_fast(photometry_fast):
    assert os.path.isfile(photometry_fast[1])


def test_photometry_mag_fast(photometry_fast):
    sn_row = photometry_fast[0][0]
    assert 15.5 < sn_row['mag'] < 15.9
    assert 0.11 < sn_row['mag_error'] < 0.18


# ------------------ EVEN SLOWER TESTS ------------------------------------------------
@pytest.mark.slow
@pytest.fixture(scope='session', autouse=True)
def photometry_slow_optical(fid=1446):
    result_table, tdir = run_phot(fid)
    yield result_table, tdir.photometry_path
    clean_up_phot(tdir.photometry_path)


@pytest.mark.slow
def test_photometry_exists_slow(photometry_slow_optical):
    assert os.path.isfile(photometry_slow_optical[1])


@pytest.mark.slow
def test_photometry_mag_slow(photometry_slow_optical):
    sn_row = photometry_slow_optical[0][0]
    assert 17.35 < sn_row['mag'] < 17.5
    assert 0.003 < sn_row['mag_error'] < 0.02


if __name__ == '__main__':
    pytest.main([__file__])
