#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Pytest fixture to create temp copy of input data, shared across all tests.

.. codeauthor:: Rasmus Handberg <rasmush@phys.au.dk>
"""

import pytest
import sys
import os
import configparser
import subprocess
import shlex

if sys.path[0] != os.path.abspath(os.path.join(os.path.dirname(__file__), '..')):
    sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))


# --------------------------------------------------------------------------------------------------
def capture_cli(script, params=None, mpiexec=False):
    if params == None:
        params = []
    if isinstance(params, str):
        params = shlex.split(params)

    cmd = [sys.executable, script.strip()] + list(params)
    if mpiexec:
        cmd = ['mpiexec', '-n', '2'] + cmd

    print(f"Command: {cmd}")
    proc = subprocess.Popen(cmd, cwd=os.path.join(os.path.dirname(__file__), '..'), stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE, universal_newlines=True)
    out, err = proc.communicate()
    exitcode = proc.returncode
    proc.kill()

    print(f"ExitCode: {exitcode:d}")
    print("StdOut:\n%s" % out.strip())
    print("StdErr:\n%s" % err.strip())
    return out, err, exitcode


# --------------------------------------------------------------------------------------------------
@pytest.fixture(scope='session')
def SETUP_CONFIG():
    """
    Fixture which sets up a dummy config-file which allows for simple testing only.
    """
    config_file = os.path.join(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'flows', 'config.ini')))
    if os.path.exists(config_file):
        yield config_file
    else:
        confstr = os.environ.get('FLOWS_CONFIG')
        if confstr is None:
            raise RuntimeError("Config file can not be set up.")

        # Write minimal config file that can be used for testing:
        config = configparser.ConfigParser()
        config.read_string(confstr)
        with open(config_file, 'w') as fid:
            config.write(fid)
            fid.flush()

        yield config_file
        os.remove(config_file)


def pytest_addoption(parser):
    parser.addoption(
        "--runslow", action="store_true", default=False, help="run slow tests"
    )


def pytest_configure(config):
    config.addinivalue_line("markers", "slow: mark test as slow to run")


def pytest_collection_modifyitems(config, items):
    if config.getoption("--runslow"):
        # --runslow given in cli: do not skip slow tests
        return
    skip_slow = pytest.mark.skip(reason="need --runslow option to run")
    for item in items:
        if "slow" in item.keywords:
            item.add_marker(skip_slow)
