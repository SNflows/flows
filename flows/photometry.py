#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Flows photometry code.

.. codeauthor:: Rasmus Handberg <rasmush@phys.au.dk>
"""

import os
import numpy as np
from bottleneck import nansum, nanmedian, allnan
from timeit import default_timer
import logging
import warnings
from copy import deepcopy

from astropy.utils.exceptions import AstropyDeprecationWarning
import astropy.units as u
import astropy.coordinates as coords
from astropy.stats import sigma_clip, SigmaClip, gaussian_fwhm_to_sigma
from astropy.table import Table, vstack
from astropy.nddata import NDData
from astropy.modeling import models, fitting
from astropy.wcs.utils import proj_plane_pixel_area

warnings.simplefilter('ignore', category=AstropyDeprecationWarning)
from photutils import DAOStarFinder, CircularAperture, CircularAnnulus, aperture_photometry
from photutils.psf import EPSFBuilder, EPSFFitter, BasicPSFPhotometry, DAOGroup, extract_stars
from photutils import Background2D, SExtractorBackground
from photutils.utils import calc_total_error

from scipy.interpolate import UnivariateSpline

from . import api
from .config import load_config
from .plots import plt, plot_image
from .version import get_version
from .load_image import load_image
from .run_imagematch import run_imagematch
from .zeropoint import bootstrap_outlier, sigma_from_Chauvenet

from .wcs import WCS

__version__ = get_version(pep440=False)

warnings.simplefilter('ignore', category=AstropyDeprecationWarning)

#--------------------------------------------------------------------------------------------------
def photometry(fileid, output_folder=None, attempt_imagematch=True):
	"""
	Run photometry.

	Parameters:
		fileid (int): File ID to process.
		output_folder (str, optional): Path to directory where output should be placed.
		attempt_imagematch (bool, optional): If no subtracted image is available, but a
			template image is, should we attempt to run ImageMatch using standard settings.
			Default=True.

	.. codeauthor:: Rasmus Handberg <rasmush@phys.au.dk>
	"""

	# Settings:
	#ref_mag_limit = 22 # Lower limit on reference target brightness
	ref_target_dist_limit = 10 * u.arcsec # Reference star must be further than this away to be included

	logger = logging.getLogger(__name__)
	tic = default_timer()

	# Use local copy of archive if configured to do so:
	config = load_config()

	# Get datafile dict from API:
	datafile = api.get_datafile(fileid)
	logger.debug("Datafile: %s", datafile)
	targetid = datafile['targetid']
	target_name = datafile['target_name']
	photfilter = datafile['photfilter']

	archive_local = config.get('photometry', 'archive_local', fallback=None)
	if archive_local is not None:
		datafile['archive_path'] = archive_local
	if not os.path.isdir(datafile['archive_path']):
		raise FileNotFoundError("ARCHIVE is not available: " + datafile['archive_path'])

	# Get the catalog containing the target and reference stars:
	# TODO: Include proper-motion to the time of observation
	catalog = api.get_catalog(targetid, output='table')
	target = catalog['target'][0]
	target_coord = coords.SkyCoord(ra=target['ra'], dec=target['decl'], unit='deg', frame='icrs')

	# Folder to save output:
	if output_folder is None:
		output_folder_root = config.get('photometry', 'output', fallback='.')
		output_folder = os.path.join(output_folder_root, target_name, '%05d' % fileid)
	logger.info("Placing output in '%s'", output_folder)
	os.makedirs(output_folder, exist_ok=True)

	# The paths to the science image:
	filepath = os.path.join(datafile['archive_path'], datafile['path'])

	# TODO: Download datafile using API to local drive:
	# TODO: Is this a security concern?
	#if archive_local:
	#	api.download_datafile(datafile, archive_local)

	# Translate photometric filter into table column:
	ref_filter = {
		'up': 'u_mag',
		'gp': 'g_mag',
		'rp': 'r_mag',
		'ip': 'i_mag',
		'zp': 'z_mag',
		'B': 'B_mag',
		'V': 'V_mag',
		'J': 'J_mag',
		'H': 'H_mag',
		'K': 'K_mag',
	}.get(photfilter, None)

	if ref_filter is None:
		logger.warning("Could not find filter '%s' in catalogs. Using default gp filter.", photfilter)
		ref_filter = 'g_mag'

	# Check that there actually are reference stars in that filter:
	if allnan(catalog['references'][ref_filter]):
		raise ValueError("No reference stars found in current photfilter.")

	# Load the image from the FITS file:
	image = load_image(filepath)

	references = extract_references(catalog, image.wcs, ref_filter, hsize, ref_target_dist_limit)

	# Calculate the targets position in the image:
	target_pixel_pos = image.wcs.all_world2pix([[target['ra'], target['decl']]], 0)[0]

	#==============================================================================================
	# BARYCENTRIC CORRECTION OF TIME
	#==============================================================================================

	ltt_bary = image.obstime.light_travel_time(target_coord, ephemeris='jpl')
	image.obstime = image.obstime.tdb + ltt_bary

	#==============================================================================================
	# BACKGROUND ESITMATION
	#==============================================================================================

	fig, ax = plt.subplots(1, 2, figsize=(20, 18))
	plot_image(image.clean, ax=ax[0], scale='log', cbar='right', title='Image')
	plot_image(image.mask, ax=ax[1], scale='linear', cbar='right', title='Mask')
	fig.savefig(os.path.join(output_folder, 'original.png'), bbox_inches='tight')
	plt.close(fig)

	# Estimate image background:
	# Not using image.clean here, since we are redefining the mask anyway
	bkg = Background2D(image.clean, (128, 128), filter_size=(5, 5),
		sigma_clip=SigmaClip(sigma=3.0),
		bkg_estimator=SExtractorBackground(),
		exclude_percentile=50.0)
	image.background = bkg.background
	image.std = bkg.background_rms_median

	# Create background-subtracted image:
	image.subclean = image.clean - image.background

	# Plot background estimation:
	fig, ax = plt.subplots(1, 3, figsize=(20, 6))
	plot_image(image.clean, ax=ax[0], scale='log', title='Original')
	plot_image(image.background, ax=ax[1], scale='log', title='Background')
	plot_image(image.subclean, ax=ax[2], scale='log', title='Background subtracted')
	fig.savefig(os.path.join(output_folder, 'background.png'), bbox_inches='tight')
	plt.close(fig)

	# TODO: Is this correct?!
	image.error = calc_total_error(image.clean, bkg.background_rms, 1.0)

	#==============================================================================================
	# DETECTION OF STARS AND MATCHING WITH CATALOG
	#==============================================================================================

	logger.info("References:\n%s", references)

	radius = 10
	fwhm_guess = 6.0
	fwhm_min = 3.5
	fwhm_max = 18.0

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
		xmin = max(x - radius, 0)
		xmax = min(x + radius + 1, image.shape[1])
		ymin = max(y - radius, 0)
		ymax = min(y + radius + 1, image.shape[0])

		curr_star = deepcopy(image.subclean[ymin:ymax, xmin:xmax])

		edge = np.zeros_like(curr_star, dtype='bool')
		edge[(0,-1),:] = True
		edge[:,(0,-1)] = True
		curr_star -= nanmedian(curr_star[edge])
		curr_star /= np.max(curr_star)

		ypos, xpos = np.mgrid[:curr_star.shape[0], :curr_star.shape[1]]
		gfit = gfitter(g2d, x=xpos, y=ypos, z=curr_star)

		fwhms[i] = gfit.x_fwhm

	masked_fwhms = np.ma.masked_array(fwhms, ~np.isfinite(fwhms))
	fwhm = np.mean(sigma_clip(masked_fwhms, maxiters=20, sigma=2.0))
	logger.info("FWHM: %f", fwhm)
	if np.isnan(fwhm):
		raise Exception("Could not estimate FWHM")

	# Use DAOStarFinder to search the image for stars, and only use reference-stars where a
	# star was actually detected close to the references-star coordinate:
	min_references = 6
	cleanout_references = (len(references) > 6)
	logger.debug("Number of references before cleaning: %d", len(references))
	if cleanout_references:
		# Arguments for DAOStarFind, starting with the strictest and ending with the
		# least strict settings to try:
		# We will stop with the first set that yield more than the minimum number
		# of reference stars.
		daofind_args = [{
			'threshold': 7 * image.std,
			'fwhm': fwhm,
			'exclude_border': True,
			'sharphi': 0.8,
			'sigma_radius': 1.1,
			'peakmax': image.peakmax
		}, {
			'threshold': 3 * image.std,
			'fwhm': fwhm,
			'roundlo': -0.5,
			'roundhi': 0.5
		}]

		# Loop through argument sets for DAOStarFind:
		for kwargs in daofind_args:

			# Run DAOStarFind with the given arguments:
			daofind_tbl = DAOStarFinder(**kwargs).find_stars(image.subclean, mask=image.mask)

			def match_catalog_to_image_with_wcs(wcs):

				references = extract_references(catalog, wcs, ref_filter, hsize, ref_target_dist_limit)

				xy_stars_grid = np.tile(daofind_tbl['xcentroid'], (len(references), 1)), np.tile(daofind_tbl['ycentroid'], (len(references), 1))
				xy_references = np.array([references['pixel_column'], references['pixel_row']]).reshape([2, len(references), 1])
				distances = np.sqrt(np.power(xy_stars_grid - xy_references, 2).sum(axis=0))
				idx_refs, idx_stars = np.any(distances <= fwhm/4, axis=1)

				return references, (idx_refs, idx_stars)

			# indices of matching reference stars and image stars
			idx_refs, idx_stars = match_catalog_to_image_with_wcs(image.wcs)[1]

			while True:

				wcs_new = WCS.from_astropy_wcs(image.wcs)
				wcs_new.adjust_with_points(
						xy = list(zip(daofind_tbl['xcentroid'][idx_stars], daofind_tbl['ycentroid'][idx_stars])),
						rd = list(zip(references['ra'][idx_refs], references['decl'][idx_refs]))
				)
				wcs_new = wcs_new.astropy_wcs

				references_new, (idx_refs_new, idx_stars_new) = match_catalog_to_image_with_wcs(wcs_new)

				# break out if wcs_new did not locate more stars
				if len(idx_refs_new) <= len(idx_refs):
					break

				references, image.wcs = references_new, wcs_new
				idx_refs, idx_stars = idx_refs_new, idx_stars_new

			indx_good = idx_refs

			logger.debug("Number of references after cleaning: %d", len(indx_good))

			if len(indx_good) >= min_references:

				references = references[indx_good]
				target_pixel_pos = image.wcs.all_world2pix([[target['ra'], target['decl']]], 0)[0]

				break

	logger.debug("Number of references after cleaning: %d", len(references))

	# Create plot of target and reference star positions:
	fig, ax = plt.subplots(1, 1, figsize=(20, 18))
	plot_image(image.subclean, ax=ax, scale='log', cbar='right', title=target_name)
	ax.scatter(references['pixel_column'], references['pixel_row'], c='r', marker='o', alpha=0.6)
	if cleanout_references:
		ax.scatter(daofind_tbl['xcentroid'], daofind_tbl['ycentroid'], c='g', marker='o', alpha=0.6)
	ax.scatter(target_pixel_pos[0], target_pixel_pos[1], marker='+', s=20, c='r')
	fig.savefig(os.path.join(output_folder, 'positions.png'), bbox_inches='tight')
	plt.close(fig)

	#==============================================================================================
	# CREATE EFFECTIVE PSF MODEL
	#==============================================================================================

	# Make cutouts of stars using extract_stars:
	# Scales with FWHM
	size = int(np.round(29*fwhm/6))
	if size % 2 == 0:
		size += 1 # Make sure it's a uneven number
	size = max(size, 15) # Never go below 15 pixels
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
				plot_image(stars[imgnr], ax=ax[i], scale='log', cmap='viridis')
			imgnr += 1

		fig.savefig(os.path.join(output_folder, 'epsf_stars%02d.png' % (k+1)), bbox_inches='tight')
		plt.close(fig)

	# Build the ePSF:
	epsf = EPSFBuilder(
		oversampling=1.0,
		maxiters=500,
		fitter=EPSFFitter(fit_boxsize=2*fwhm),
		progress_bar=True
	)(stars)[0]

	logger.info('Successfully built PSF model')

	fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 15))
	plot_image(epsf.data, ax=ax1, cmap='viridis')

	fwhms = []
	bad_epsf_detected = False
	for a, ax in ((0, ax3), (1, ax2)):
		# Collapse the PDF along this axis:
		profile = epsf.data.sum(axis=a)
		itop = profile.argmax()
		poffset = profile[itop]/2

		# Run a spline through the points, but subtract half of the peak value, and find the roots:
		# We have to use a cubic spline, since roots() is not supported for other splines
		# for some reason
		profile_intp = UnivariateSpline(np.arange(0, len(profile)), profile - poffset, k=3, s=0, ext=3)
		lr = profile_intp.roots()

		# Plot the profile and spline:
		x_fine = np.linspace(-0.5, len(profile)-0.5, 500)
		ax.plot(profile, 'k.-')
		ax.plot(x_fine, profile_intp(x_fine) + poffset, 'g-')
		ax.axvline(itop)
		ax.set_xlim(-0.5, len(profile)-0.5)

		# Do some sanity checks on the ePSF:
		# It should pass 50% exactly twice and have the maximum inside that region.
		# I.e. it should be a single gaussian-like peak
		if len(lr) != 2 or itop < lr[0] or itop > lr[1]:
			logger.error("Bad PSF along axis %d", a)
			bad_epsf_detected = True
		else:
			axis_fwhm = lr[1] - lr[0]
			fwhms.append(axis_fwhm)
			ax.axvspan(lr[0], lr[1], facecolor='g', alpha=0.2)

	# Save the ePSF figure:
	ax4.axis('off')
	fig.savefig(os.path.join(output_folder, 'epsf.png'), bbox_inches='tight')
	plt.close(fig)

	# There was a problem with the ePSF:
	if bad_epsf_detected:
		raise Exception("Bad ePSF detected.")

	# Let's make the final FWHM the largest one we found:
	fwhm = np.max(fwhms)
	logger.info("Final FWHM based on ePSF: %f", fwhm)

	#==============================================================================================
	# COORDINATES TO DO PHOTOMETRY AT
	#==============================================================================================

	coordinates = np.array([[ref['pixel_column'], ref['pixel_row']] for ref in references])

	# Add the main target position as the first entry for doing photometry directly in the
	# science image:
	coordinates = np.concatenate(([target_pixel_pos], coordinates), axis=0)

	#==============================================================================================
	# APERTURE PHOTOMETRY
	#==============================================================================================

	# Define apertures for aperture photometry:
	apertures = CircularAperture(coordinates, r=fwhm)
	annuli = CircularAnnulus(coordinates, r_in=1.5*fwhm, r_out=2.5*fwhm)

	apphot_tbl = aperture_photometry(image.subclean, [apertures, annuli], mask=image.mask, error=image.error)

	logger.debug("Aperture Photometry Table:\n%s", apphot_tbl)
	logger.info('Apperature Photometry Success')

	#==============================================================================================
	# PSF PHOTOMETRY
	#==============================================================================================

	# Are we fixing the postions?
	epsf.fixed.update({'x_0': False, 'y_0': False})

	# Create photometry object:
	photometry = BasicPSFPhotometry(
		group_maker=DAOGroup(fwhm),
		bkg_estimator=SExtractorBackground(),
		psf_model=epsf,
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

	#==============================================================================================
	# TEMPLATE SUBTRACTION AND TARGET PHOTOMETRY
	#==============================================================================================

	# Find the pixel-scale of the science image:
	pixel_area = proj_plane_pixel_area(image.wcs.celestial)
	pixel_scale = np.sqrt(pixel_area)*3600 # arcsec/pixel
	#print(image.wcs.celestial.cunit) % Doesn't work?
	logger.info("Science image pixel scale: %f", pixel_scale)

	diffimage = None
	if datafile.get('diffimg') is not None:

		diffimg_path = os.path.join(datafile['archive_path'], datafile['diffimg']['path'])
		diffimage = load_image(diffimg_path)
		diffimage = diffimage.image

	elif attempt_imagematch and datafile.get('template') is not None:
		# Run the template subtraction, and get back
		# the science image where the template has been subtracted:
		diffimage = run_imagematch(datafile, target, star_coord=coordinates, fwhm=fwhm, pixel_scale=pixel_scale)

	# We have a diff image, so let's do photometry of the target using this:
	if diffimage is not None:
		# Include mask from original image:
		diffimage = np.ma.masked_array(diffimage, image.mask)

		# Create apertures around the target:
		apertures = CircularAperture(target_pixel_pos, r=fwhm)
		annuli = CircularAnnulus(target_pixel_pos, r_in=1.5*fwhm, r_out=2.5*fwhm)

		# Create two plots of the difference image:
		fig, ax = plt.subplots(1, 1, squeeze=True, figsize=(20, 20))
		plot_image(diffimage, ax=ax, cbar='right', title=target_name)
		ax.plot(target_pixel_pos[0], target_pixel_pos[1], marker='+', color='r')
		fig.savefig(os.path.join(output_folder, 'diffimg.png'), bbox_inches='tight')
		apertures.plot(color='r')
		annuli.plot(color='k')
		ax.set_xlim(target_pixel_pos[0]-50, target_pixel_pos[0]+50)
		ax.set_ylim(target_pixel_pos[1]-50, target_pixel_pos[1]+50)
		fig.savefig(os.path.join(output_folder, 'diffimg_zoom.png'), bbox_inches='tight')
		plt.close(fig)

		# Run aperture photometry on subtracted image:
		target_apphot_tbl = aperture_photometry(diffimage, [apertures, annuli], mask=image.mask, error=image.error)

		# Run PSF photometry on template subtracted image:
		target_psfphot_tbl = photometry(
			diffimage,
			init_guesses=Table(target_pixel_pos, names=['x_0', 'y_0'])
		)

		# Combine the output tables from the target and the reference stars into one:
		apphot_tbl = vstack([target_apphot_tbl, apphot_tbl], join_type='exact')
		psfphot_tbl = vstack([target_psfphot_tbl, psfphot_tbl], join_type='exact')

	# Build results table:
	tab = references.copy()
	tab.insert_row(0, {'starid': 0, 'ra': target['ra'], 'decl': target['decl'], 'pixel_column': target_pixel_pos[0], 'pixel_row': target_pixel_pos[1]})
	if diffimage is not None:
		tab.insert_row(0, {'starid': -1, 'ra': target['ra'], 'decl': target['decl'], 'pixel_column': target_pixel_pos[0], 'pixel_row': target_pixel_pos[1]})
	indx_main_target = (tab['starid'] <= 0)
	for key in ('pm_ra', 'pm_dec', 'gaia_mag', 'gaia_bp_mag', 'gaia_rp_mag', 'B_mag', 'V_mag', 'H_mag','J_mag','K_mag', 'u_mag', 'g_mag', 'r_mag', 'i_mag', 'z_mag'):
		for i in np.where(indx_main_target)[0]: # No idea why this is needed, but giving a boolean array as slice doesn't work
			tab[i][key] = np.NaN

	# Subtract background estimated from annuli:
	flux_aperture = apphot_tbl['aperture_sum_0'] - (apphot_tbl['aperture_sum_1'] / annuli.area()) * apertures.area()
	flux_aperture_error = np.sqrt(apphot_tbl['aperture_sum_err_0']**2 + (apphot_tbl['aperture_sum_err_1']/annuli.area() * apertures.area())**2)

	# Add table columns with results:
	tab['flux_aperture'] = flux_aperture/image.exptime
	tab['flux_aperture_error'] = flux_aperture_error/image.exptime
	tab['flux_psf'] = psfphot_tbl['flux_fit']/image.exptime
	tab['flux_psf_error'] = psfphot_tbl['flux_unc']/image.exptime
	tab['pixel_column_psf_fit'] = psfphot_tbl['x_fit']
	tab['pixel_row_psf_fit'] = psfphot_tbl['y_fit']
	tab['pixel_column_psf_fit_error'] = psfphot_tbl['x_0_unc']
	tab['pixel_row_psf_fit_error'] = psfphot_tbl['y_0_unc']

	# Check that we got valid photometry:
	if np.any(~np.isfinite(tab[indx_main_target]['flux_psf'])) or np.any(~np.isfinite(tab[indx_main_target]['flux_psf_error'])):
		raise Exception("Target magnitude is undefined.")

	#==============================================================================================
	# CALIBRATE
	#==============================================================================================

	# Convert PSF fluxes to magnitudes:
	mag_inst = -2.5 * np.log10(tab['flux_psf'])
	mag_inst_err = (2.5/np.log(10)) * (tab['flux_psf_error'] / tab['flux_psf'])

	# Corresponding magnitudes in catalog:
	#TODO: add color terms here
	mag_catalog = tab[ref_filter]

	# Mask out things that should not be used in calibration:
	use_for_calibration = np.ones_like(mag_catalog, dtype='bool')
	use_for_calibration[indx_main_target] = False # Do not use target for calibration
	use_for_calibration[~np.isfinite(mag_inst) | ~np.isfinite(mag_catalog)] = False

	# Just creating some short-hands:
	x = mag_catalog[use_for_calibration]
	y = mag_inst[use_for_calibration]
	yerr = mag_inst_err[use_for_calibration]
	weights = 1.0/yerr**2

	# Fit linear function with fixed slope, using sigma-clipping:
	model = models.Linear1D(slope=1, fixed={'slope': True})
	fitter = fitting.FittingWithOutlierRemoval(fitting.LinearLSQFitter(), sigma_clip, sigma=3.0)
	best_fit, sigma_clipped = fitter(model, x, y, weights=weights)

	# Extract zero-point and estimate its error using a single weighted fit:
	# I don't know why there is not an error-estimate attached directly to the Parameter?
	zp = -1*best_fit.intercept.value # Negative, because that is the way zeropoints are usually defined

	weights[sigma_clipped] = 0 # Trick to make following expression simpler
	N = len(weights.nonzero()[0])
	if N > 1:
		zp_error = np.sqrt( N * nansum(weights*(y - best_fit(x))**2) / nansum(weights) / (N-1) )
	else:
		zp_error = np.NaN
	logger.info('Leastsquare ZP = {0:0.3f}, ZP_error = {1:0.3f}'.format(zp, zp_error))

	# Determine sigma clipping sigma according to Chauvenet method
	# But don't allow less than sigma = sigmamin, setting to 1.5 for now.
	# Should maybe be 2?
	sigmamin = 1.5
	sigChauv = sigma_from_Chauvenet(len(x))
	sigChauv = sigChauv if sigChauv >= sigmamin else sigmamin

	# Extract zero point and error using bootstrap method
	Nboot = 1000
	logger.info('Running bootstrap with sigma = {0:0.2f} and n = {1:0.0f}'.format(sigChauv,Nboot))
	pars = bootstrap_outlier(x, y, yerr, n=Nboot, model=model, fitter=fitting.LinearLSQFitter,
		outlier=sigma_clip, outlier_kwargs={'sigma':sigChauv}, summary='median',
		error='bootstrap', return_vals=False)

	zp_bs = pars['intercept'] * -1.0
	zp_error_bs = pars['intercept_error']

	logger.info('Bootstrapped ZP = {0:0.3f}, ZP_error = {1:0.3f}'.format(zp_bs,zp_error_bs))

	# Check that difference is not large
	zp_diff = 0.4
	if np.abs(zp_bs - zp) >= zp_diff:
		logger.warning("Bootstrap and weighted LSQ ZPs differ by {:0.2f}, \
		which is more than the allowed {:0.2f} mag.".format(np.abs(zp_bs - zp), zp_diff))

	# Add calibrated magnitudes to the photometry table:
	tab['mag'] = mag_inst + zp_bs
	tab['mag_error'] = np.sqrt(mag_inst_err**2 + zp_error_bs**2)

	fig, ax = plt.subplots(1, 1)
	ax.errorbar(x, y, yerr=yerr, fmt='k.')
	ax.scatter(x[sigma_clipped], y[sigma_clipped], marker='x', c='r')
	ax.plot(x, best_fit(x), color='g', linewidth=3)
	ax.set_xlabel('Catalog magnitude')
	ax.set_ylabel('Instrumental magnitude')
	fig.savefig(os.path.join(output_folder, 'calibration.png'), bbox_inches='tight')
	plt.close(fig)

	# Check that we got valid photometry:
	if not np.isfinite(tab[0]['mag']) or not np.isfinite(tab[0]['mag_error']):
		raise Exception("Target magnitude is undefined.")

	#==============================================================================================
	# SAVE PHOTOMETRY
	#==============================================================================================

	# Descriptions of columns:
	tab['flux_aperture'].unit = u.count/u.second
	tab['flux_aperture_error'].unit = u.count/u.second
	tab['flux_psf'].unit = u.count/u.second
	tab['flux_psf_error'].unit = u.count/u.second
	tab['pixel_column'].unit = u.pixel
	tab['pixel_row'].unit = u.pixel
	tab['pixel_column_psf_fit'].unit = u.pixel
	tab['pixel_row_psf_fit'].unit = u.pixel
	tab['pixel_column_psf_fit_error'].unit = u.pixel
	tab['pixel_row_psf_fit_error'].unit = u.pixel

	# Meta-data:
	tab.meta['version'] = __version__
	tab.meta['fileid'] = fileid
	tab.meta['template'] = None if datafile.get('template') is None else datafile['template']['fileid']
	tab.meta['diffimg'] = None if datafile.get('diffimg') is None else datafile['diffimg']['fileid']
	tab.meta['photfilter'] = photfilter
	tab.meta['fwhm'] = fwhm * u.pixel
	tab.meta['pixel_scale'] = pixel_scale * u.arcsec/u.pixel
	tab.meta['seeing'] = (fwhm*pixel_scale) * u.arcsec
	tab.meta['obstime-bmjd'] = float(image.obstime.mjd)
	tab.meta['zp'] = zp_bs
	tab.meta['zp_error'] = zp_error_bs
	tab.meta['zp_diff'] = np.abs(zp_bs - zp)
	tab.meta['zp_error_weights'] = zp_error

	# Filepath where to save photometry:
	photometry_output = os.path.join(output_folder, 'photometry.ecsv')

	# Write the final table to file:
	tab.write(photometry_output, format='ascii.ecsv', delimiter=',', overwrite=True)

	toc = default_timer()

	logger.info("------------------------------------------------------")
	logger.info("Success!")
	logger.info("Main target: %f +/- %f", tab[0]['mag'], tab[0]['mag_error'])
	logger.info("Photometry took: %f seconds", toc-tic)

	return photometry_output

#--------------------------------------------------------------------------------------------------

def extract_references(catalog, wcs, ref_filter='g_mag', hsize=10, ref_target_dist_limit=10*u.arcsec):

	references, target = catalog['references'], catalog['target'][0]
	references.sort(ref_filter)

	references['pixel_column'], references['pixel_row'] = x, y = \
			list(zip(*wcs.all_world2pix(list(zip(references['ra'], references['decl'])), 0)))

	targ_coord = coords.SkyCoord(ra=target['ra'], dec=target['decl'], unit='deg', frame='icrs')
	refs_coord = coords.SkyCoord(ra=references['ra'], dec=references['decl'], unit='deg', frame='icrs')
	references = references[
			(targ_coord.separation(refs_coord) > ref_target_dist_limit) &
			(x > hsize) & (x < (image.shape[1] - 1 - hsize)) &
			(y > hsize) & (y < (image.shape[0] - 1 - hsize))
	]
	# 		& (references[ref_filter] < ref_mag_limit)

	return references
