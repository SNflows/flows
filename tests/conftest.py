#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Pytest fixture to create temp copy of input data, shared across all tests.

.. codeauthor:: Rasmus Handberg <rasmush@phys.au.dk>
"""

#import pytest
import sys
import os.path
#import shutil

if sys.path[0] != os.path.abspath(os.path.join(os.path.dirname(__file__), '..')):
	sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
