#!/usr/bin/env python
# -*- coding: utf-8 -*-
# flake8: noqa

from .photometry import photometry
from .catalogs import download_catalog
from .visibility import visibility
from .config import load_config

from .version import get_version

__version__ = get_version(pep440=False)
