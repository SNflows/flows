#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Flows photometry code.

.. codeauthor:: Rasmus Handberg <rasmush@phys.au.dk>
"""

import numpy as np
import astropy.units as u
import astropy.coordinates as coords
from astropy.io import fits
from astropy.time import Time
from astropy.wcs import WCS, FITSFixedWarning
from . import api
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
	with fits.open(FILENAME, mode='readonly', memmap=True) as hdul:
		hdr = hdul[0].header
		origin = hdr.get('ORIGIN')

		image.image = np.asarray(hdul[0].data)
		image.shape = image.image.shape

		if origin == 'LCOGT':
			image.mask = np.asarray(hdul['BPM'].data, dtype='bool')
		else:
			image.mask = np.zeros_like(image.image, dtype='bool')

		image.mask |= ~np.isfinite(image.image)

		image.clean = np.ma.masked_array(image.image, image.mask)

		# World Coordinate System:
		with warnings.catch_warnings():
			warnings.simplefilter('ignore', category=FITSFixedWarning)
			image.wcs = WCS(hdr)

		# Specific headers everyone can agree on:
		image.exptime = float(hdr['EXPTIME']) # * u.second
		image.peakmax = None # Maximum value above which data is not to be trusted

		# Timestamp:
		if origin == 'LCOGT':
			sites = api.sites.get_all_sites()
			site_keywords = {s['site_keyword']: s for s in sites}
			image.site = site_keywords.get(hdr['SITE'], None)

			observatory = coords.EarthLocation.from_geodetic(lat=hdr['LATITUDE'], lon=hdr['LONGITUD'], height=hdr['HEIGHT'])
			image.obstime = Time(hdr['MJD-OBS'], format='mjd', scale='utc', location=observatory)

			image.photfilter = hdr['FILTER']

			# Get non-linear limit
			# TODO: Use actual or some fraction of the non-linearity limit
			#image.peakmax = hdr.get('MAXLIN') # Presumed non-linearity limit from header
			image.peakmax = 60000 #From experience, this one is better.

		elif hdr.get('TELESCOP') == 'NOT' and hdr.get('INSTRUME') in ('ALFOSC FASU', 'ALFOSC_FASU') and hdr.get('OBS_MODE') == 'IMAGING':
			image.site = api.get_site(5) # Hard-coded the siteid for NOT
			image.obstime = Time(hdr['DATE-AVG'], format='isot', scale='utc', location=image.site['EarthLocation'])

			image.photfilter = {
				'B Bes': 'B',
				'V Bes': 'V',
				'g SDSS': 'gp',
				'r SDSS': 'rp',
				'i SDSS': 'ip',
				'u SDSS': 'up'
			}.get(hdr['FILTER'].replace('_', ' '), hdr['FILTER'])

			# Get non-linear limit
			# Obtained from http://www.not.iac.es/instruments/detectors/CCD14/LED-linearity/20181026-200-1x1.pdf
			# TODO: grab these from a table for all detector setups of ALFOSC
			image.peakmax = 80000 # For ALFOSC D, 1x1, 200; the standard for SNe.

		elif hdr.get('FPA.TELESCOPE') == 'PS1' and hdr.get('FPA.INSTRUMENT') == 'GPC1':
			image.site = api.get_site(6) # Hard-coded the siteid for Pan-STARRS1
			image.obstime = Time(hdr['MJD-OBS'], format='mjd', scale='utc', location=image.site['EarthLocation'])

			image.photfilter = {
				'g.00000': 'gp',
				'r.00000': 'rp',
				'i.00000': 'ip',
				'z.00000': 'zp'
			}.get(hdr['FPA.FILTER'], hdr['FPA.FILTER'])

		elif hdr.get('TELESCOP') == 'Liverpool Telescope':
			# Liverpool telescope
			image.site = api.get_site(8) # Hard-coded the siteid for Liverpool Telescope
			image.obstime = Time(hdr['DATE-OBS'], format='isot', scale='utc', location=image.site['EarthLocation'])
			image.obstime += 0.5*image.exptime * u.second # Make time centre of exposure
			image.photfilter = {
				'Bessel-B': 'B',
				'Bessel-V': 'V',
				'SDSS-U': 'up',
				'SDSS-G': 'gp',
				'SDSS-R': 'rp',
				'SDSS-I': 'ip',
				'SDSS-Z': 'zp'
			}.get(hdr['FILTER1'], hdr['FILTER1'])

		else:
			raise Exception("Could not determine origin of image")

	return image
