#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Flows photometry code.

.. codeauthor:: Rasmus Handberg <rasmush@phys.au.dk>
"""

import os
import numpy as np
from bottleneck import nanstd
from timeit import default_timer
import logging
import sys

import astropy.units as u
import astropy.coordinates as coords
from astropy.time import Time
from astropy.wcs import WCS
from astropy.io import fits
from astropy.stats import sigma_clip, gaussian_sigma_to_fwhm, SigmaClip, gaussian_fwhm_to_sigma
from astropy.table import Table
from astropy.nddata import NDData
from astropy.modeling import models, fitting

from photutils import DAOStarFinder, CircularAperture, CircularAnnulus, aperture_photometry
from photutils.psf import EPSFBuilder, EPSFFitter, BasicPSFPhotometry, DAOGroup, extract_stars
from photutils import Background2D, SExtractorBackground

from scipy.interpolate import UnivariateSpline

from .catalogs import get_catalog
from .datafiles import get_datafile
from .plots import plt, plot_image

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

		image.image = hdul[0].data
		image.wcs = WCS(hdr)
		image.mask = np.asarray(hdul[2].data, dtype='bool')
		image.clean = np.ma.masked_array(image.image, image.mask)
		image.shape = image.image.shape

		# Specific headers:
		image.exptime = float(hdul[0].header['EXPTIME']) # * u.second

		observatory = coords.EarthLocation.from_geodetic(lat=hdr['LATITUDE'], lon=hdr['LONGITUD'], height=hdr['HEIGHT'])
		image.obstime = Time(hdr['MJD-OBS'], format='mjd', scale='utc', location=observatory)

	return image

#--------------------------------------------------------------------------------------------------
def photometry(fileid=None):
	"""


	.. codeauthor:: Rasmus Handberg <rasmush@phys.au.dk>
	"""

	plt.switch_backend('Qt5Agg')

	logger = logging.getLogger(__name__)
	tic = default_timer()

	# Get datafile dict from API:
	datafile = get_datafile(fileid)
	print(datafile)
	FILENAME = os.path.join(r'C:\Users\au195407\Documents\flows_archive', datafile['path']) # datafile['archive_path']
	targetid = datafile['targetid']
	photfilter = datafile['photfilter']

	# Settings:
	background_cutoff = 1000 # All pixels above this threshold are masked during background estimation
	ref_mag_limit = 17 # Lower limit on reference target brightness
	ref_target_dist_limit = 30 # Reference star must be further than this away to be included

	# Translate photometric filter into table column:
	if photfilter in ('B', 'V', 'gp'):
		ref_filter = 'g_mag'
	elif photfilter == 'rp':
		ref_filter = 'r_mag'
	elif photfilter == 'ip':
		ref_filter = 'i_mag'
	elif photfilter == 'zp':
		ref_filter = 'z_mag'
	else:
		ref_filter = 'g_mag'

	# Load the image from the FITS file:
	image = load_image(FILENAME)

	# Get the catalog containing the target and reference stars:
	# TODO: Include proper-motion to the time of observation
	catalog = get_catalog(targetid, output='table')
	target = catalog['target'][0]
	references = catalog['references']
	references.sort(ref_filter)

	# Extract information about target:
	target_name = str(target['target_name'])
	target_coord = coords.SkyCoord(ra=target['ra'], dec=target['decl'], unit='deg', frame='icrs')

	# Folder to save output:
	output_folder = os.path.join(os.path.dirname(FILENAME), '%04d' % fileid)
	os.makedirs(output_folder, exist_ok=True)

	# Calculate pixel-coordinates of references:
	row_col_coords = image.wcs.all_world2pix(np.array([[ref['ra'], ref['decl']] for ref in references]), 0)
	references['pixel_column'] = row_col_coords[:,0]
	references['pixel_row'] = row_col_coords[:,1]

	# Calculate the targets position in the image:
	target_pixel_pos = image.wcs.all_world2pix([[target['ra'], target['decl']]], 0)[0]

	# Clean out the references:
	hsize = 10
	x = references['pixel_column']
	y = references['pixel_row']
	references = references[(np.sqrt((x - target_pixel_pos[0])**2 + (y - target_pixel_pos[1])**2) > ref_target_dist_limit)
		& (references[ref_filter] < ref_mag_limit)
		& (x > hsize) & (x < (image.shape[1] - 1 - hsize))
		& (y > hsize) & (y < (image.shape[0] - 1 - hsize))]

	#==============================================================================================
	# BARYCENTRIC CORRECTION OF TIME
	#==============================================================================================

	ltt_bary = image.obstime.light_travel_time(target_coord, ephemeris='jpl')
	image.obstime = image.obstime.tdb + ltt_bary

	#==============================================================================================
	# BACKGROUND ESITMATION
	#==============================================================================================

	# Estimate image background:
	# Not using image.clean here, since we are redefining the mask anyway
	bkg = Background2D(image.image, (64, 64), filter_size=(5, 5),
		mask=image.mask | (image.clean > background_cutoff),
		sigma_clip=SigmaClip(sigma=3.0),
		bkg_estimator=SExtractorBackground()
		)
	image.background = bkg.background

	# Create background-subtracted image:
	image.subclean = image.clean - image.background

	# Plot background estimation:
	fig, ax = plt.subplots(1, 3, figsize=(20, 6))
	plot_image(image.clean, ax=ax[0], scale='log', title='Original')
	plot_image(image.background, ax=ax[1], scale='log', title='Background')
	plot_image(image.subclean, ax=ax[2], scale='log', title='Background subtracted')
	fig.savefig(os.path.join(output_folder, 'background.png'))

	#==============================================================================================
	# DETECTION OF STARS AND MATCHING WITH CATALOG
	#==============================================================================================

	logger.info("References:\n%s", references)

	radius=10
	fwhm_guess = 6.0
	fwhm_min = 3.5
	fwhm_max = 13.5

	# Extract stars sub-images:
	#stars = extract_stars(
	#	NDData(data=image.subclean, mask=image.mask),
	#	stars_for_epsf,
	#	size=size
	#)

	# Set up 2D Gaussian model for fitting to reference stars:
	g2d = models.Gaussian2D(amplitude=1.0, x_mean=radius, y_mean=radius, x_stddev=fwhm_guess*gaussian_fwhm_to_sigma)
	g2d.amplitude.bounds = (0.1, 2.0)
	g2d.x_mean.bounds = (0.5*radius, 1.5*radius)
	g2d.y_mean.bounds = (0.5*radius, 1.5*radius)
	g2d.x_stddev.bounds = (fwhm_min * gaussian_fwhm_to_sigma, fwhm_max * gaussian_fwhm_to_sigma)
	g2d.y_stddev.tied = lambda model: model.x_stddev
	g2d.theta.fixed = True

	gfitter = fitting.LevMarLSQFitter()

	fwhms = np.full(len(references), np.NaN)
	for i, (x, y) in enumerate(zip(references['pixel_column'], references['pixel_row'])):
		x = int(np.round(x))
		y = int(np.round(y))
		x0, y0, width, height = x - radius, y - radius, 2 * radius, 2 * radius
		cutout = slice(y0 - 1, y0 + height), slice(x0 - 1, x0 + width)

		curr_star = image.subclean[cutout] / np.max(image.subclean[cutout])
		npix = len(curr_star)

		ypos, xpos = np.mgrid[:npix, :npix]
		gfit = gfitter(g2d, x=xpos, y=ypos, z=curr_star)

		fwhms[i] = gfit.x_fwhm

	mask = ~np.isfinite(fwhms) | (fwhms <= fwhm_min) | (fwhms >= fwhm_max)
	masked_fwhms = np.ma.MaskedArray(fwhms, mask)

	fwhm = np.mean(sigma_clip(masked_fwhms, maxiters=20, sigma=2.0))
	logger.info("FWHM: %f", fwhm)

	# Use DAOStarFinder to search the image for stars, and only use reference-stars where a
	# star was actually detected close to the references-star coordinate:
	daofind_tbl = DAOStarFinder(100, fwhm=fwhm, roundlo=-0.5, roundhi=0.5).find_stars(image.subclean, mask=image.mask)
	indx_good = np.zeros(len(references), dtype='bool')
	for k, ref in enumerate(references):
		dist = np.sqrt( (daofind_tbl['xcentroid'] - ref['pixel_column'])**2 + (daofind_tbl['ycentroid'] - ref['pixel_row'])**2 )
		if np.any(dist <= fwhm/4): # Cutoff set somewhat arbitrery
			indx_good[k] = True

	references = references[indx_good]

	fig, ax = plt.subplots(1, 1, figsize=(20, 18))
	plot_image(image.subclean, ax=ax, scale='log', make_cbar=True, title=target_name)
	ax.scatter(references['pixel_column'], references['pixel_row'], c='r', alpha=0.3)
	ax.scatter(daofind_tbl['xcentroid'], daofind_tbl['ycentroid'], c='g', alpha=0.3)
	ax.scatter(target_pixel_pos[0], target_pixel_pos[1], marker='+', c='r')
	fig.savefig(os.path.join(output_folder, 'positions.png'))
	plt.show()

	#==============================================================================================
	# CREATE EFFECTIVE PSF MODEL
	#==============================================================================================

	# Make cutouts of stars using extract_stars:
	size = 29 # TODO: Scale with fwhm
	hsize = (size - 1) / 2

	x = references['pixel_column']
	y = references['pixel_row']
	mask_near_edge = ((x > hsize) & (x < (image.shape[1] - 1 - hsize))
		& (y > hsize) & (y < (image.shape[0] - 1 - hsize)))

	stars_for_epsf = Table()
	stars_for_epsf['x'] = x[mask_near_edge]
	stars_for_epsf['y'] = y[mask_near_edge]

	# Store which stars were used in ePSF in the table:
	logger.info("Number of stars used for ePSF: %d", len(stars_for_epsf))
	references['used_for_epsf'] = mask_near_edge

	# Extract stars sub-images:
	stars = extract_stars(
		NDData(data=image.subclean, mask=image.mask),
		stars_for_epsf,
		size=size
	)

	# Plot the stars being used for ePSF:
	nrows = 5
	ncols = 5
	imgnr = 0
	for k in range(int(np.ceil(len(stars_for_epsf)/(nrows*ncols)))):
		fig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=(20, 20), squeeze=True)
		ax = ax.ravel()
		for i in range(nrows*ncols):
			if imgnr > len(stars_for_epsf)-1:
				ax[i].axis('off')
			else:
				plot_image(stars[imgnr], ax=ax[i], scale='log', cmap='viridis', xlabel=None, ylabel=None)
			imgnr += 1

		fig.savefig(os.path.join(output_folder, 'epsf_stars%02d.png' % (k+1)))
		plt.close(fig)

	# Build the ePSF:
	image.epsf = EPSFBuilder(
		oversampling=1.0,
		maxiters=500,
		fitter=EPSFFitter(fit_boxsize=2*fwhm),
		progress_bar=True
	)(stars)[0]

	logger.info('Successfully built PSF model')

	fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 15))
	plot_image(image.epsf.data, ax=ax1, cmap='viridis', xlabel=None, ylabel=None)

	fwhms = []
	for a, ax in ((0, ax3), (1, ax2)):
		# Collapse the PDF along this axis:
		profile = image.epsf.data.sum(axis=a)
		itop = profile.argmax()
		poffset = profile[itop]/2

		# Run a spline through the points, but subtract half of the peak value, and find the roots:
		# We have to use a cubic spline, since roots() is not supported for other splines
		# for some reason
		profile_intp = UnivariateSpline(np.arange(0, len(profile)), profile - poffset, k=3, s=0, ext=3)
		lr = profile_intp.roots()
		axis_fwhm = lr[1] - lr[0]

		fwhms.append(axis_fwhm)

		x_fine = np.linspace(-0.5, len(profile)-0.5, 500)

		ax.plot(profile, 'k.-')
		ax.plot(x_fine, profile_intp(x_fine) + poffset, 'g-')
		ax.axvline(itop)
		ax.axvspan(lr[0], lr[1], facecolor='g', alpha=0.2)
		ax.set_xlim(-0.5, len(profile)-0.5)

	# Let's make the final FWHM the largest one we found:
	fwhm = np.max(fwhms)
	logger.info("Final FWHM based on ePSF: %f", fwhm)

	#ax2.axvspan(itop - fwhm/2, itop + fwhm/2, facecolor='b', alpha=0.2)
	#ax3.axvspan(itop - fwhm/2, itop + fwhm/2, facecolor='b', alpha=0.2)
	ax4.axis('off')

	fig.savefig(os.path.join(output_folder, 'epsf.png'))
	plt.show()

	#==============================================================================================
	# COORDINATES TO DO PHOTOMETRY AT
	#==============================================================================================

	coordinates = np.array([[ref['pixel_column'], ref['pixel_row']] for ref in references])

	# Add the main target position as the first entry:
	coordinates = np.concatenate(([target_pixel_pos], coordinates), axis=0)

	#==============================================================================================
	# APERTURE PHOTOMETRY
	#==============================================================================================

	# Define apertures for aperture photometry:
	apertures = CircularAperture(coordinates, r=fwhm)
	annuli = CircularAnnulus(coordinates, r_in=1.5*fwhm, r_out=2.5*fwhm)

	apphot_tbl = aperture_photometry(image.subclean, [apertures, annuli], mask=image.mask)

	logger.debug("Aperture Photometry Table:\n%s", apphot_tbl)

	# Subtract background estimated from annuli:
	flux_aperture = apphot_tbl['aperture_sum_0'] - (apphot_tbl['aperture_sum_1'] / annuli.area()) * apertures.area()
	#flux_aperture_error = np.sqrt(error**2 + (phot['aperture_sum_err_1']/annuli.area() * apertures.area() )**2)
	flux_aperture_error = 0

	logger.info('Apperature Photometry Success')

	#==============================================================================================
	# PSF PHOTOMETRY
	#==============================================================================================

	# Are we fixing the postions?
	image.epsf.fixed.update({'x_0': False, 'y_0': False})

	# Create photometry object:
	photometry = BasicPSFPhotometry(
		group_maker=DAOGroup(fwhm),
		bkg_estimator=SExtractorBackground(),
		psf_model=image.epsf,
		fitter=fitting.LevMarLSQFitter(),
		fitshape=size,
		aperture_radius=fwhm
	)

	psfphot_tbl = photometry(
		image=image.subclean,
		init_guesses=Table(coordinates, names=['x_0', 'y_0'])
	)

	logger.debug("PSF Photometry Table:\n%s", psfphot_tbl)
	logger.info('PSF Photometry Success')

	# Build results table:
	tab = references.copy()
	tab.insert_row(0, {'starid': 0, 'ra': target['ra'], 'decl': target['decl'], 'pixel_column': target_pixel_pos[0], 'pixel_row': target_pixel_pos[1]})
	tab[0]['H_mag'] = None
	#

	tab['flux_aperture'] = flux_aperture
	tab['flux_aperture_error'] = flux_aperture_error
	tab['flux_psf'] = psfphot_tbl['flux_fit']
	tab['flux_psf_error'] = psfphot_tbl['flux_unc']
	tab['pixel_column_psf_fit'] = psfphot_tbl['x_fit']
	tab['pixel_row_psf_fit'] = psfphot_tbl['y_fit']
	tab['pixel_column_psf_fit_error'] = psfphot_tbl['x_0_unc']
	tab['pixel_row_psf_fit_error'] = psfphot_tbl['y_0_unc']

	#==============================================================================================
	# CALIBRATE
	#==============================================================================================

	mag_catalog = tab[ref_filter]
	mag_inst = -2.5 * np.log10(tab['flux_psf'] / image.exptime)
	mag_inst_err = tab['flux_psf_error'] / tab['flux_psf']/image.exptime * 1.0857 # ??????????????
	print(mag_inst_err)

	# Mask out things that should not be used in calibration:
	use_for_calibration = np.ones_like(mag_catalog, dtype='bool')
	use_for_calibration[0] = False # Do not use target for calibration
	use_for_calibration[~np.isfinite(mag_inst) | ~np.isfinite(mag_catalog)] = False

	x = mag_catalog[use_for_calibration]
	y = mag_inst[use_for_calibration]
	yerr = mag_inst_err[use_for_calibration]

	# Fit linear function with fixed slope, using sigma-clipping:
	model = models.Linear1D(slope=1, fixed={'slope': True})
	fitter = fitting.FittingWithOutlierRemoval(fitting.LinearLSQFitter(), sigma_clip, sigma=3.0)
	best_fit, sigma_clipped = fitter(model, x, y, weights=1.0/yerr**2)

	# Extract zero-point and estimate its error:
	# I don't know why there is not an error-estimate attached directly to the Parameter?
	zp = best_fit.intercept.value
	zp_error = nanstd(y[~sigma_clipped] - best_fit(x[~sigma_clipped]))

	# Add calibrated magnitudes to the photometry table:
	tab['mag'] = mag_inst - zp
	tab['mag_error'] = np.sqrt(mag_inst_err**2 + zp_error**2)

	fig, ax = plt.subplots(1, 1)
	ax.errorbar(x, y, yerr=yerr, fmt='k.')
	ax.scatter(x[sigma_clipped], y[sigma_clipped], marker='x', c='r')
	ax.plot(x, best_fit(x), color='g', linewidth=3)
	ax.set_xlabel('Catalog magnitude')
	ax.set_ylabel('Instrumental magnitude')
	fig.savefig(os.path.join(output_folder, 'calibration.png'))
	plt.show()

	#==============================================================================================
	# SAVE PHOTOMETRY
	#==============================================================================================

	# Descriptions of columns:
	tab['pixel_column'].unit = u.pixel
	tab['pixel_row'].unit = u.pixel
	tab['pixel_column_psf_fit'].unit = u.pixel
	tab['pixel_row_psf_fit'].unit = u.pixel
	tab['pixel_column_psf_fit_error'].unit = u.pixel
	tab['pixel_row_psf_fit_error'].unit = u.pixel

	# Meta-data:
	tab.meta['fileid'] = fileid
	tab.meta['photfilter'] = photfilter
	tab.meta['fwhm'] = fwhm
	tab.meta['obstime-bmjd'] = float(image.obstime.mjd)
	tab.meta['zp'] = zp
	tab.meta['zp_error'] = zp_error

	print(tab)
	tab.write(os.path.join(output_folder, 'photometry.ecsv'), format='ascii.ecsv', delimiter=',', overwrite=True)

	toc = default_timer()
	logger.info("Photometry took: %f seconds", toc-tic)
