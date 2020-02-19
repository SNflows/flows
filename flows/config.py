#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

"""

import os.path
import configparser

#--------------------------------------------------------------------------------------------------
def load_config():

	config_file = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'config.ini')

	config = configparser.ConfigParser()
	config.read(config_file)
	return config