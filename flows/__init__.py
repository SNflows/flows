"""
FLOWS pipeline package
"""
from .photometry import do_phot as photometry
from .catalogs import download_catalog
from .visibility import visibility
from .version import get_version

__version__ = get_version(pep440=False)
