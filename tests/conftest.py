#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Pytest fixture to create temp copy of input data, shared across all tests.

.. codeauthor:: Rasmus Handberg <rasmush@phys.au.dk>
"""

import pytest
import sys
import os
#import shutil
import configparser

if sys.path[0] != os.path.abspath(os.path.join(os.path.dirname(__file__), '..')):
	sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

#--------------------------------------------------------------------------------------------------
@pytest.fixture(scope='session')
def SETUP_CONFIG():
	"""
	Fixture which sets up a dummy config-file which allows for simple testing only.
	"""
	config_file = os.path.join(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'flows', 'config.ini')))
	if os.path.exists(config_file):
		yield config_file
	else:
		token = os.environ.get('FLOWS_CONFIG_TOKEN')
		if token is None:
			raise RuntimeError("Config file can not be set up.")

		# Write minimal config file that can be used for testing:
		config = configparser.ConfigParser()
		config['api'] = {'token': token}
		with open(config_file, 'w') as fid:
			config.write(fid)
			fid.flush()

		yield config_file
		os.remove(config_file)
