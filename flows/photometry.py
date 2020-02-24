#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import numpy as np
import matplotlib.pyplot as plt

from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from astropy.io import fits
from astropy.stats import sigma_clip, sigma_clipped_stats, gaussian_sigma_to_fwhm
from astropy.table import Table
from astropy.nddata import NDData
from astropy.modeling.fitting import LevMarLSQFitter

from photutils import DAOStarFinder, CircularAperture, CircularAnnulus, aperture_photometry, EPSFBuilder, EPSFFitter
from photutils.psf import BasicPSFPhotometry, DAOGroup, extract_stars
from photutils.background import MMMBackground

from imexam.imexamine import Imexamine
#import pandas as pd
#import calviacat as cvc

from .catalogs import get_catalog
from .aadc_db import AADC_DB
from .plots import plt, plot_image

#--------------------------------------------------------------------------------------------------
def load_image(FILENAME):

	# Get image and WCS, find stars, remove galaxies
	image = type('image', (object,), dict()) # image container

	# get image and wcs solution
	with fits.open(FILENAME, mode='readonly', memmap=True) as hdul:
		image.image = hdul[0].data
		image.wcs = WCS(hdul[0].header)
		image.mask = np.asarray(hdul[2].data, dtype='bool')
		image.clean = np.ma.masked_array(image.image, image.mask)
		image.shape = image.image.shape

	return image

#--------------------------------------------------------------------------------------------------
def photometry(fileid=None, target=2, FWHM=6):

	plt.switch_backend('Qt5Agg')

	if fileid is not None:
		with AADC_DB() as db:
			db.cursor.execute("SELECT files.*, files_archives.path AS archive_path FROM flows.files INNER JOIN files_archives ON files.archive=files_archives.archive WHERE fileid=%s;", [fileid])
			row = db.cursor.fetchone()

			FILENAME = os.path.join(r'C:\Users\au195407\Documents\flows_archive', row['path']) # row['archive_path']
			target = row['targetid']

	print(FILENAME)
	print(target)

	image = load_image(FILENAME)

	catalog = get_catalog(target)

	# Extract information about target:
	#target_name = catalog['target']['target_name']
	#ra = catalog['target']['ra']
	#dec = catalog['target']['decl']
	#coord_SN = SkyCoord(ra, dec, unit=('deg', 'deg'), frame='icrs')

	# References:
	references = Table(rows=catalog['references'])
	references.sort('J_mag')

	row_col_coords = image.wcs.all_world2pix(np.array([[ref['ra'], ref['decl']] for ref in references]), 0)
	references['pixel_column'] = row_col_coords[:,0]
	references['pixel_row'] = row_col_coords[:,1]

	# Clean out the references:
	hsize = 10
	ref_mag_limit = 15
	x = references['pixel_column']
	y = references['pixel_row']
	references = references[(references['J_mag'] < ref_mag_limit)
		& (x > hsize) & (x < (image.shape[1] - 1 - hsize))
		& (y > hsize) & (y < (image.shape[0] - 1 - hsize))]

	# get find sources in the image
	image.background, std = sigma_clipped_stats(image.image, sigma=3, maxiters=10, mask=image.mask)[1:]
	image.subclean = image.clean - image.background

	print(references)

	plt.figure()
	plot_image(image.subclean, scale='log', make_cbar=True)
	plt.scatter(references['pixel_column'], references['pixel_row'], c='r', alpha=0.3)
	plt.show()

	FWHM = 6 # pixels, First Guess, will be derived in script.
	plots = Imexamine()
	fwhms = np.full(len(references), np.NaN)
	for i, (x, y) in enumerate(zip(references['pixel_column'], references['pixel_row'])):
		try:
			fwhms[i] = (plots.line_fit(x, y, image.subclean, genplot=False).stddev_0 * gaussian_sigma_to_fwhm)
		except ValueError:
			fwhms[i] = np.NaN

	mask = ~np.isfinite(fwhms) | (fwhms < 3.5) | (fwhms > 11.0)
	masked_fwhms = np.ma.MaskedArray(fwhms, mask)

	fwhms_clean = sigma_clip(masked_fwhms, maxiters=20, sigma=2.0)
	fwhm = int(np.round(np.mean(fwhms_clean)))
	print(fwhm)

	"""
	X, Y, Z = *np.meshgrid(*map(np.arange, image.image.shape[::-1])), image.image - image.background
	for x, y in zip(np.round(image.x).astype(int), np.round(image.y).astype(int)):
		i = slice(y - FWHM, y + FWHM + 1), slice(x - FWHM, x + FWHM + 1)
		#x, y = np.round( (Z[i] * (X[i], Y[i])).sum(axis=(1,2)) / Z[i].sum() ).astype(int)
		#i = slice(y - FWHM, y + FWHM + 1), slice(x - FWHM, x + FWHM + 1)
		image.stars.append(i)

	local_sequence=np.where(~fwhms_clean.mask)[0]


	if saveseq:
		np.savetxt(NEW_FILENAME, list(zip(image.ra[local_sequence], image.dec[local_sequence])))


	# ra, dec

	FWHM = fwhm # pixels
	ANNULUS_RADII = FWHM+2, FWHM+5 # pixels
	OVERSAMPLING = 1 # factor


	coordinates = type('coordinates', (object,), dict())
	coordinates.ra, coordinates.dec = np.loadtxt(COORDINATES, unpack=True, usecols=(0,1))
	coordinates.x, coordinates.y = zip(*image.wcs.wcs_world2pix(list(zip(*(coordinates.ra, coordinates.dec))), 0))


	image.stars = [
		(slice(y - FWHM, y + FWHM + 1), slice(x - FWHM, x + FWHM + 1))
		for x, y in zip(np.round(coordinates.x).astype(int), np.round(coordinates.y).astype(int))
	]

	image.centered_stars = image.stars.copy()
	coordinates.centered_x = list(coordinates.x).copy()
	coordinates.centered_y = list(coordinates.y).copy()



	apertures = CircularAperture(zip(coordinates.centered_x, coordinates.centered_y), FWHM)
	annuli = CircularAnnulus(zip(coordinates.centered_x, coordinates.centered_y), *ANNULUS_RADII)

	#fig = plt.figure(figsize=(20, 2 * (len(image.stars) // 10 + 1)))
	#for i in range(len(image.stars)):

		#aperture, annulus = apertures.to_mask()[i], annuli.to_mask()[i]

		#ax = fig.add_subplot(len(image.stars) // 10 + 1, 10, i + 1)
		#ax.text(0, 0, i, color='black', fontsize=14, transform=ax.transAxes)

		#img = aperture.cutout(image.image)
		#img = plt.get_cmap('viridis')((img - img.min()) / (img.max() - img.min()))
		#img[:,:,3] = aperture.data
		#plt.imshow(img, origin='lower', extent=aperture.bbox.extent)

		#img = annulus.cutout(image.image)
		#img = plt.get_cmap('viridis')((img - img.min()) / (img.max() - img.min()))
		#img[:,:,3] = annulus.data
		#plt.imshow(img, origin='lower', extent=annulus.bbox.extent)

		#apertures.plot(indices=i, color='white')

	#fig.tight_layout()
	"""

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
	print("Number of stars used for ePSF: %d" % len(stars_for_epsf))
	references['used_for_epsf'] = mask_near_edge

	print(stars_for_epsf)

	# Extract stars sub-images:
	stars = extract_stars(
		NDData(data=image.subclean, mask=image.mask),
		stars_for_epsf,
		size=size
	)

	# Plot the stars being used for ePSF:
	nrows = 5
	ncols = 5
	fig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=(20, 20), squeeze=True)
	ax = ax.ravel()
	for i in range(nrows*ncols):
		plot_image(stars[i], ax=ax[i], scale='log', cmap='viridis', xlabel=None, ylabel=None)

	plt.show(block=True)

	# Build the ePSF:
	image.epsf = EPSFBuilder(
		oversampling=4.0,
		maxiters=500,
		fitter=EPSFFitter(fit_boxsize=2*fwhm),
		progress_bar=True
	)(stars)[0]

	plt.figure()
	plot_image(image.epsf.data, cmap='viridis', make_cbar=True)
	plt.show()

	print('Successfully built PSF model')

	profile = image.epsf.data.sum(axis=0)
	itop = profile.argmax()

	left = np.argmin((profile[:itop] - profile.max()/2)**2)
	right = np.argmin((profile[itop:] - profile.max()/2)**2) + itop

	FWHM=np.ceil((right - left) / OVERSAMPLING).astype(int)

	FWHM=int(FWHM)

	#APERTURE PHOT
	coordinates.aperture = list()

	error = image.background + image.image - image.background
	photometry = aperture_photometry(image.image, [apertures, annuli], error=np.sqrt(image.image))

	for phot in photometry:

		flux, error = phot['aperture_sum_0'], phot['aperture_sum_err_0']

		flux = flux - phot['aperture_sum_1'] / annuli.area() * apertures.area()
		error = np.sqrt(error**2 + (phot['aperture_sum_err_1'] * apertures.area() / annuli.area())**2)

		coordinates.aperture.append((flux, error))

	print('App. Phot Success')
	#PSF PHOT

	coordinates.psf = list()

	# Create photometry object:
	photometry = BasicPSFPhotometry(
		group_maker = DAOGroup(FWHM),
		bkg_estimator = MMMBackground(),
		psf_model = image.epsf,
		fitter = LevMarLSQFitter(),
		fitshape = 29,
		aperture_radius = FWHM
	)

	for i, star in enumerate(image.centered_stars):

		image.epsf.fixed.update(dict(x_0=False, y_0=False))

		res = photometry(
			image = image.image - image.background,
			init_guesses=Table([
				dict(x_0=coordinates.centered_x[i], y_0=coordinates.centered_y[i])
			])
		)
		coordinates.psf.append((res['flux_fit'].data[0], res['flux_unc'].data[0]))

	print('Psf Phot Success')

	data = np.vstack([coordinates.ra, coordinates.dec, *zip(*coordinates.psf), *zip(*coordinates.aperture)])
	#header = 'aperture = %d, annulus = %d, %d' % (FWHM, *ANNULUS_RADII)
	header ='ra dec psf psf_err ap ap_err'
	np.savetxt(PHOT_FILENAME, np.transpose(data), header=header,comments='#')

	#==============================================================================================
	# CALIBRATE
	#==============================================================================================
	"""
	df=pd.read_csv(PHOT_FILENAME,delim_whitespace=True,escapechar='#')
	AT=Table.from_pandas(df)
	c=SkyCoord(AT['ra'],AT['dec'],unit=('deg', 'deg'), frame='icrs')
	i_rm=coord_SN.match_to_catalog_sky(c)[0] #Location of SN

	wsid='1844987713'
	password='jeknun-dosxaw-7cAvdo'

	ps1 = cvc.PanSTARRS1(path+'cat/cat'+fname_alt+'.db')
	ps1.fetch_field(c)

	objids, distances = ps1.xmatch(c)


	i_inst = -2.5 * np.log10(AT['psf']/image.header['EXPTIME'])
	i_err = AT['psf_err'] / AT['psf']/image.header['EXPTIME'] * 1.0857

	filtername=filter_dictionary[image.header['FILTER']]

	zp_mean,zp_med,zp_err,m,gmi_constant=ps1.cal_constant(objids, i_inst,filtername,mlim=[13,18])

	F_inst=(i_inst+zp_med)
	F_inst_err=np.sqrt(i_err**2 + zp_err**2)
	"""

	#==============================================================================================
	# SAVE PHOTOMETRY
	#==============================================================================================
	"""
	mjd=image.header['MJD-OBS']
	row_dict={'mjd':mjd,'mag':F_inst[i_rm],'err':F_inst_err[i_rm],'filter':filtername}

	final_phot_file=path+'phot/photdf.csv'
	if os.path.exists(final_phot_file):
	    photdf=pd.read_csv(final_phot_file)

	    if (~np.any(photdf.mjd==mjd)):
	    #if True:
	        photdf_newrow=pd.DataFrame({'mjd':mjd,'mag':F_inst[i_rm],'err':F_inst_err[i_rm],'filter':filtername},index=[0])
	        photdf=photdf.append(photdf_newrow)
	    else:
	        row_change=np.where(photdf.mjd==mjd)[0][0]
	        photdf_newrow=pd.DataFrame({'mjd':mjd,'mag':F_inst[i_rm],'err':F_inst_err[i_rm],'filter':filtername},index=[row_change])

	        photdf.iloc[row_change]=(mjd,F_inst[i_rm],F_inst_err[i_rm],filtername)

	else:
	    photdf=pd.DataFrame({'mjd':mjd,'mag':F_inst[i_rm],'err':F_inst_err[i_rm],'filter':filtername},index=[0])

	photdf.to_csv(final_phot_file,index=False)
	photdf=pd.read_csv(final_phot_file)
	photdf.sort_values('mjd')
	"""
