#!/usr/bin/env python
# -*- coding: utf-8 -*-

from .photometry import photometry
from .download_catalog import download_catalog

from .version import get_version
__version__ = get_version(pep440=False)
