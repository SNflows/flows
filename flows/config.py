#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

.. codeauthor:: Rasmus Handberg <rasmush@phys.au.dk>
"""

import os.path
import configparser
from functools import lru_cache

#--------------------------------------------------------------------------------------------------
@lru_cache
def load_config():
	"""
	Load configuration file.

	Returns:
		``configparser.ConfigParser``: Configuration file.

	.. codeauthor:: Rasmus Handberg <rasmush@phys.au.dk>
	"""

	config_file = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'config.ini')
	if not os.path.isfile(config_file):
		raise FileNotFoundError("config.ini file not found")

	config = configparser.ConfigParser()
	config.read(config_file)
	return config
