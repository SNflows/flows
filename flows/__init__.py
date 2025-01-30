"""
FLOWS pipeline package
"""
from .photometry import timed_photometry as photometry
from .catalogs import download_catalog, delete_catalog
from .visibility import visibility
from .version import get_version
from .load_image import load_image
from .instruments import INSTRUMENTS

__version__ = get_version(pep440=False)
