#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Flows photometry code.

.. codeauthor:: Rasmus Handberg <rasmush@phys.au.dk>
"""

import numpy as np
import astropy.coordinates as coords
from astropy.io import fits
from astropy.time import Time
from astropy.wcs import WCS, FITSFixedWarning
import warnings

#--------------------------------------------------------------------------------------------------
def load_image(FILENAME):
	"""
	Load FITS image.

	Parameters:
		FILENAME (string): Path to FITS file to be loaded.

	Returns:
		object: Image constainer.

	.. codeauthor:: Rasmus Handberg <rasmush@phys.au.dk>
	"""

	# Get image and WCS, find stars, remove galaxies
	image = type('image', (object,), dict()) # image container

	# get image and wcs solution
	with fits.open(FILENAME, mode='readonly') as hdul:
		hdr = hdul[0].header

		image.image = np.asarray(hdul[0].data)
		image.mask = np.asarray(hdul[2].data, dtype='bool')
		image.clean = np.ma.masked_array(image.image, image.mask)
		image.shape = image.image.shape

		# World Coordinate System:
		with warnings.catch_warnings():
			warnings.simplefilter('ignore', category=FITSFixedWarning)
			image.wcs = WCS(hdr)

		# Specific headers:
		image.exptime = float(hdr['EXPTIME']) # * u.second

		observatory = coords.EarthLocation.from_geodetic(lat=hdr['LATITUDE'], lon=hdr['LONGITUD'], height=hdr['HEIGHT'])
		image.obstime = Time(hdr['MJD-OBS'], format='mjd', scale='utc', location=observatory)

	return image
