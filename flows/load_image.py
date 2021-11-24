#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Flows photometry code.

.. codeauthor:: Emir Karamehmetoglu <emir.k@phys.au.dk>
.. codeauthor:: Rasmus Handberg <rasmush@phys.au.dk>
"""

import numpy as np
import warnings
import logging
import astropy.units as u
import astropy.coordinates as coords
from astropy.io import fits
from astropy.time import Time
from astropy.wcs import WCS, FITSFixedWarning
from flows import api
from dataclasses import dataclass, field
import typing
from abc import ABC, abstractmethod


# --------------------------------------------------------------------------------------------------
@dataclass
class FlowsImage:
	image: np.ndarray
	header: typing.Dict
	mask: np.ndarray = None
	peakmax: float = None

	def __post_init__(self):
		self.shape = self.image.shape
		self.wcs = self.create_wcs()
		# Make empty mask
		if not self.mask:
			self.mask = np.zeros_like(self.image, dtype='bool')
		self.check_finite()

	def check_finite(self):
		self.mask |= ~np.isfinite(self.image)

	def update_mask(self, mask):
		self.mask = mask
		self.check_finite()

	# def get_observatory_params(self, instrument=None):
	# 		if not instrument:
	# 		instrument = self.instrument
	# 	site = instrument.get_site()
	# 	obstime = instrument.get_obstime()
	# 	photfilter = instrument.get_photfilter()

	def create_wcs(self) -> WCS:
		with warnings.catch_warnings():
			warnings.simplefilter('ignore', category=FITSFixedWarning)
			return WCS(header=self.header, relax=True)

	def create_masked_image(self):
		"""Warning: this is destructive and will overwrite image data setting masked values to NaN"""
		self.image[self.mask] = np.NaN
		self.clean = np.ma.masked_array(data=self.image, mask=self.mask, copy=False)


class AbstractInstrument(ABC):
	peakmax: int = None

	@abstractmethod
	def __init__(self):
		pass

	@abstractmethod
	def get_site(self):
		pass

	@abstractmethod
	def get_exptime(self):
		pass

	@abstractmethod
	def get_obstime(self):
		pass

	@abstractmethod
	def get_photfilter(self):
		pass


class Instrument(AbstractInstrument):
	peakmax: int = None
	siteid: int = None

	def __init__(self, image: FlowsImage):
		self.image = image
		self.image.peakmax = self.peakmax
		self.image.site = self.get_site()
		self.image.exptime = self.get_exptime()
		self.image.obstime = self.get_obstime()
		self.image.photfilter = self.get_photfilter()
		self.image.create_masked_image()

	def get_site(self):
		if self.siteid is not None:
			return api.get_site(self.siteid)

	def get_exptime(self):
		exptime = self.image.header.get('EXPTIME', None)
		if exptime is None:
			raise ValueError("Image exposure time could not be extracted")
		return exptime

	def get_obstime(self):
		'''Default for JD, jd, utc.'''
		return Time(self.image.header['JD'], format='jd', scale='utc',
					location=self.image.site['EarthLocation'])

	def get_photfilter(self):
		return self.image.header['FILTER']


class LCOGT(Instrument):
	peakmax: int = 60000

	def get_site(self):
		sites = api.sites.get_all_sites()
		site_keywords = {s['site_keyword']: s for s in sites}
		site = site_keywords.get(self.image.header['SITE'], None)
		return site

	def get_obstime(self):
		observatory = coords.EarthLocation.from_geodetic(lat=self.image.header['LATITUDE'],
														 lon=self.image.header['LONGITUD'],
														 height=self.image.header['HEIGHT'])
		obstime = Time(self.image.header['DATE-OBS'], format='isot', scale='utc', location=observatory)
		obstime += 0.5 * self.image.exptime * u.second  # Make time centre of exposure
		return obstime

	def get_photfilter(self):
		photfilter = {'zs': 'zp'}.get(self.image.header['FILTER'], self.image.header['FILTER'])
		return photfilter


class HAWKI(Instrument):
	siteid = 2  # Hard-coded the siteid for ESO Paranal, VLT, UT4

	def __init__(self, image: FlowsImage):
		super().__init__(image)
		self.get_obtype()

	def get_obstime(self):
		obstime = Time(self.image.header['DATE-OBS'], format='isot', scale='utc',
								  location=self.image.site['EarthLocation'])
		obstime += 0.5 * self.image.exptime * u.second  # Make time centre of exposure
		return obstime

	def get_obtype(self):
		ob_type = self.image.header["HIERARCH ESO OCS DET1 IMGNAME"].split('_')[-1]
		if "Auto" in ob_type:
			self.image.ob_type = 'Autojitter'
		elif "Fixed" in ob_type:
			self.image.ob_type = 'FixedOffset'
		else:
			raise RuntimeError("Image OB Type not AutoJitter or FixedOffset")


def get_image_extension(hdul: fits.HDUList, target_coord: coords.SkyCoord = None, fallback_extension: int = None):
	# For HAWKI multi-extension images we search the extensions for which one contains
	# the target, Create Image from that extension.
	target_radec = [[target_coord.icrs.ra.deg, target_coord.icrs.dec.deg]]

	for k in range(1, 5):
		w = WCS(header=hdul[k].header, relax=True)
		s = [hdul[k].header['NAXIS2'], hdul[k].header['NAXIS1']]
		pix = w.all_world2pix(target_radec, 0).flatten()
		if -0.5 <= pix[0] <= s[1] - 0.5 and -0.5 <= pix[1] <= s[0] - 0.5:
			return k
	if not fallback_extension is None:
		return fallback_extension
	else:
		raise RuntimeError(f"Could not find image extension that target is on!")


class ALFOSC(Instrument):
	# Obtained from http://www.not.iac.es/instruments/detectors/CCD14/LED-linearity/20181026-200-1x1.pdf
	peakmax = 80000  # For ALFOSC D, 1x1, 200; the standard for SNe.
	siteid = 5

	def get_obstime(self):
		return Time(self.image.header['DATE-AVG'], format='isot', scale='utc',
					location=self.image.site['EarthLocation'])

	def get_photfilter(self):
		# Sometimes data from NOT does not have the FILTER keyword,
		# in which case we have to try to figure out which filter
		# was used based on some of the other headers:
		if 'FILTER' in self.image.header:
			photfilter = {
				'B Bes': 'B',
				'V Bes': 'V',
				'R Bes': 'R',
				'g SDSS': 'gp',
				'r SDSS': 'rp',
				'i SDSS': 'ip',
				'i int': 'ip',  # Interference filter
				'u SDSS': 'up',
				'z SDSS': 'zp'
			}.get(self.image.header['FILTER'].replace('_', ' '), self.image.header['FILTER'])
		else:
			filters_used = []
			for check_headers in ('ALFLTNM', 'FAFLTNM', 'FBFLTNM'):
				isopen = self.image.header.get(check_headers).strip().lower() != 'open'
				if self.image.header.get(check_headers) and isopen:
					filters_used.append(self.image.header.get(check_headers).strip())
			if len(filters_used) == 1:
				photfilter = {
					'V_Bes 530_80': 'V',
					'R_Bes 650_130': 'R',
					"g'_SDSS 480_145": 'gp',
					"r'_SDSS 618_148": 'rp',
					"i'_SDSS 771_171": 'ip',
					'i_int 797_157': 'ip',  # Interference filter
					"z'_SDSS 832_LP": 'zp'
				}.get(filters_used[0].replace('  ', ' '), filters_used[0])
			else:
				raise RuntimeError("Could not determine filter used.")

		return photfilter


class NOTCAM(Instrument):
	siteid = 5

	def get_obstime(self):
		return Time(self.image.header['DATE-AVG'], format='isot', scale='utc',
					location=self.image.site['EarthLocation'])

	def get_photfilter(self):
		# Does NOTCAM data sometimes contain a FILTER header?
		# if not we have to try to figure out which filter
		# was used based on some of the other headers:
		if 'FILTER' in self.image.header:
			raise RuntimeError("NOTCAM: Filter keyword defined")
		filters_used = []
		for check_headers in ('NCFLTNM1', 'NCFLTNM2'):
			isopen = self.image.header.get(check_headers).strip().lower() != 'open'
			if self.image.header.get(check_headers) and isopen:
				filters_used.append(self.image.header.get(check_headers).strip())
		if len(filters_used) == 1:
			photfilter = {}.get(filters_used[0], filters_used[0])
		else:
			raise RuntimeError("Could not determine filter used.")
		return photfilter


class PS1(Instrument):
	siteid = 6

	def get_obstime(self):
		return Time(self.image.header['MJD-OBS'], format='mjd', scale='utc',
					location=self.image.site['EarthLocation'])

	def get_photfilter(self):
		photfilter = {
			'g.00000': 'gp',
			'r.00000': 'rp',
			'i.00000': 'ip',
			'z.00000': 'zp'
		}.get(self.image.header['FPA.FILTER'], self.image.header['FPA.FILTER'])
		return photfilter


class Liverpool(Instrument):
	siteid = 8

	def get_obstime(self):
		obstime = Time(self.image.header['DATE-OBS'], format='isot', scale='utc',
					   location=self.image.site['EarthLocation'])
		obstime += 0.5 * self.image.exptime * u.second  # Make time centre of exposure
		return obstime

	def get_photfilter(self):
		photfilter = {
			'Bessel-B': 'B',
			'Bessell-B': 'B',
			'Bessel-V': 'V',
			'Bessell-V': 'V',
			'SDSS-U': 'up',
			'SDSS-G': 'gp',
			'SDSS-R': 'rp',
			'SDSS-I': 'ip',
			'SDSS-Z': 'zp'
		}.get(self.image.header['FILTER1'], self.image.header['FILTER1'])
		return photfilter


class Omega2000(Instrument):
	siteid = 9

	def get_obstime(self):
		obstime = Time(self.image.header['MJD-OBS'], format='mjd', scale='utc',
					   location=self.image.site['EarthLocation'])
		obstime += 0.5 * self.image.exptime * u.second
		return obstime


class Swope(Instrument):
	siteid = 10

	def get_photfilter(self):
		photfilter = {
			'u': 'up',
			'g': 'gp',
			'r': 'rp',
			'i': 'ip',
		}.get(self.image.header['FILTER'], self.image.header['FILTER'])
		return photfilter


class Dupont(Instrument):
	siteid = 14

	def get_photfilter(self):
		photfilter = {
			'u': 'up',
			'g': 'gp',
			'r': 'rp',
			'i': 'ip',
		}.get(self.image.header['FILTER'], self.image.header['FILTER'])
		return photfilter


class RetroCam(Instrument):
	siteid = 16

	def get_photfilter(self):
		photfilter = {
			'Yc': 'Y',
			'Hc': 'H',
			'Jo': 'J',
		}.get(self.image.header['FILTER'], self.image.header['FILTER'])
		return photfilter


class Baade(Instrument):
	siteid = 11

	def get_exptime(self):
		exptime = super().get_exptime()
		exptime *= int(self.image.header['NCOMBINE'])  # EXPTIME is only for a single exposure
		return exptime

	def get_photfilter(self):
		photfilter = {
			'Ks': 'K',
			'J1': 'Y',
		}.get(self.image.header['FILTER'], self.image.header['FILTER'])
		return photfilter


class Sofi(Instrument):
	siteid = 12

	def get_obstime(self):
		if 'TMID' in self.image.header:
			obstime = Time(self.image.header['TMID'], format='mjd', scale='utc',
						   location=self.image.site['EarthLocation'])
		else:
			obstime = Time(self.image.header['MJD-OBS'], format='mjd', scale='utc',
						   location=self.image.site['EarthLocation'])
			obstime += 0.5 * self.image.exptime * u.second  # Make time centre of exposure
		return obstime

	def get_photfilter(self):
		hdr = self.image.header
		photfilter_translate = {
			'Ks': 'K'
		}
		if 'FILTER' in hdr:
			photfilter = photfilter_translate.get(hdr['FILTER'], hdr['FILTER'])
		else:
			filters_used = []
			for check_headers in ('ESO INS FILT1 ID', 'ESO INS FILT2 ID'):
				if hdr.get(check_headers) and hdr.get(check_headers).strip().lower() != 'open':
					filters_used.append(hdr.get(check_headers).strip())
			if len(filters_used) == 1:
				photfilter = photfilter_translate.get(filters_used[0], filters_used[0])
			else:
				raise RuntimeError("Could not determine filter used.")
		return photfilter


class EFOSC(Instrument):
	siteid = 15

	def get_obstime(self):
		obstime = Time(self.image.header['DATE-OBS'], format='isot', scale='utc',
					   location=self.image.site['EarthLocation'])
		obstime += 0.5 * self.image.exptime * u.second  # Make time centre of exposure
		return obstime

	def get_photfilter(self):
		hdr = self.image.header
		photfilter = {
			'g782': 'gp',
			'r784': 'rp',
			'i705': 'ip',
			'B639': 'B',
			'V641': 'V'
		}.get(hdr['FILTER'], hdr['FILTER'])
		return photfilter


class AstroNIRCam(Instrument):
	siteid = 13

	def get_exptime(self):
		return self.image.header.get('FULL_EXP', super().get_exptime())

	def get_obstime(self):
		hdr = self.image.header
		if 'MIDPOINT' in hdr:
			obstime = Time(hdr['MIDPOINT'], format='isot', scale='utc', location=self.image.site['EarthLocation'])
		else:
			obstime = Time(hdr['MJD-AVG'], format='mjd', scale='utc', location=self.image.site['EarthLocation'])
		return obstime

	def get_photfilter(self):
		hdr = self.image.header
		photfilter = {
			'H_Open': 'H',
			'K_Open': 'K',
		}.get(hdr['FILTER'], hdr['FILTER'])
		return photfilter


class OmegaCam(Instrument):
	siteid = 18  # Hard-coded the siteid for ESO VLT Survey telescope

	def get_obstime(self):
		obstime = Time(self.image.header['MJD-OBS'], format='mjd', scale='utc',
					   location=self.image.site['EarthLocation'])
		obstime += 0.5 * self.image.exptime * u.second  # Make time centre of exposure
		return obstime

	def get_photfilter(self):
		hdr = self.image.header
		photfilter = {
			'i_SDSS': 'ip'
		}.get(hdr['ESO INS FILT1 NAME'], hdr['ESO INS FILT1 NAME'])
		return photfilter


class AndiCam(Instrument):
	siteid = 20  # Hard-coded the siteid for ANDICAM at Cerro Tololo Interamerican Observatory (CTIO)

	def get_obstime(self):
		obstime = super().get_obstime()
		obstime += 0.5 * self.image.exptime * u.second
		return obstime

	def get_photfilter(self):
		return self.image.header['CCDFLTID']


class PairTel(Instrument):
	siteid = 21

	def get_obstime(self):
		hdr = self.image.header
		time_start = Time(hdr['STRT_CPU'], format='iso', scale='utc', location=self.image.site['EarthLocation'])
		time_stop = Time(hdr['STOP_CPU'], format='iso', scale='utc', location=self.image.site['EarthLocation'])
		obstime = time_start + 0.5 * (time_stop - time_start)
		return obstime

	def get_photfilter(self):
		hdr = self.image.header
		photfilter = {
			'j': 'J',
			'h': 'H',
			'k': 'K',
		}.get(hdr['FILTER'], hdr['FILTER'])
		return photfilter

class TJO(Instrument):
	siteid = 22

	def get_obstime(self):
		obstime = super().get_obstime()
		obstime += 0.5 * self.image.exptime * u.second
		return obstime


instruments = {
	'LCOGT': LCOGT,
	'HAWKI': HAWKI,
	'ALFOSC': ALFOSC,
	'NOTCAM': NOTCAM,
	'PS1': PS1,
	'Liverpool': Liverpool,
	'Omega2000': Omega2000,
	'Swope': Swope,
	'Dupont': Dupont,
	'Retrocam': RetroCam,
	'Baade': Baade,
	'Sofi': Sofi,
	'EFOSC': EFOSC,
	'AstroNIRCam': AstroNIRCam,
	'OmegaCam': OmegaCam,
	'AndiCam': AndiCam,
	'PairTel': PairTel,
	'TJO': TJO
}


def edge_mask(img, value=0):
	"""
	Create boolean mask of given value near edge of image.

	Parameters:
		img (ndarray): Image of
		value (float): Value to detect near edge. Default=0.

	Returns:
		ndarray: Pixel mask with given values on the edge of image.

	.. codeauthor:: Rasmus Handberg <rasmush@phys.au.dk>
	"""

	mask1 = (img == value)
	mask = np.zeros_like(img, dtype='bool')

	# Mask entire rows and columns which are only the value:
	mask[np.all(mask1, axis=1), :] = True
	mask[:, np.all(mask1, axis=0)] = True

	# Detect "uneven" edges column-wise in image:
	a = np.argmin(mask1, axis=0)
	b = np.argmin(np.flipud(mask1), axis=0)
	for col in range(img.shape[1]):
		if mask1[0, col]:
			mask[:a[col], col] = True
		if mask1[-1, col]:
			mask[-b[col]:, col] = True

	# Detect "uneven" edges row-wise in image:
	a = np.argmin(mask1, axis=1)
	b = np.argmin(np.fliplr(mask1), axis=1)
	for row in range(img.shape[0]):
		if mask1[row, 0]:
			mask[row, :a[row]] = True
		if mask1[row, -1]:
			mask[row, -b[row]:] = True

	return mask


def new_load_image(filename: str, target_coord: coords.SkyCoord = None):
	"""
	Load FITS image using FlowsImage class and Instrument Classes.

	Parameters:
		filename (str): Path to FITS file to be loaded.
		target_coord (:class:`astropy.coordinates.SkyCoord`): Coordinates of target.
			Only used for HAWKI images to determine which image extension to load,
			for all other images it is ignored.

	Returns:
		FlowsImage: instance of FlowsImage with valuues populated based on instrument.

	.. codeauthor:: Emir Karamehmetoglu <emir.k@phys.au.dk>
	.. codeauthor:: Rasmus Handberg <rasmush@phys.au.dk>
	"""
	logger = logging.getLogger(__name__)
	ext = 0  # Default extension in HDUList, individual instruments may override this.
	mask = None  # Instrument can override, default is to only mask all non-finite values, override is additive.

	# Read fits image, Structural Pattern Match to specific instrument.
	with fits.open(filename, mode='readonly') as hdul:
		hdr = hdul[ext].header
		origin = hdr.get('ORIGIN', '')
		telescope = hdr.get('TELESCOP', '')
		instrument = hdr.get('INSTRUME', '')

		instrument_name = None
		# Pattern matching begins here, ideally we use 3.10 pattern matching, or a dictionary lookup.
		while instrument_name is None:
			# LCOGT
			if origin == "LCOGT":
				instrument_name = 'LCOGT'
				if 'BPM' in hdul:
					mask = np.asarray(hdul['BPM'].data, dtype='bool')
				else:
					logger.warning('LCOGT image does not contain bad pixel map. Not applying mask.')
				break
			# HAWKI
			elif origin == 'ESO-PARANAL' and telescope == 'ESO-VLT-U4' and instrument == 'HAWKI' and hdr.get(
				'PRODCATG') == 'SCIENCE.MEFIMAGE':
				instrument_name = 'HAWKI'

				if target_coord is None:
					raise ValueError("TARGET_COORD is needed for HAWKI images to find the correct extension")
				if not isinstance(target_coord, coords.SkyCoord):
					if len(target_coord) == 2:
						target_coord = coords.SkyCoord(ra=target_coord[0] * u.deg, dec=target_coord[1] * u.deg,
													   frame='icrs')
					else:
						raise ValueError("TARGET_COORD is needed for HAWKI images to find the correct extension")
				ext = get_image_extension(hdul, target_coord)  # Find the one with the SN in it.
				hdr = hdul[ext].header + hdul[0].header
				break

			# NOT - ALFOSC
			elif telescope == 'NOT' and instrument in ('ALFOSC FASU','ALFOSC_FASU') and hdr.get('OBS_MODE','').lower() == 'imaging':
				instrument_name = 'ALFOSC'
				break

			elif telescope == 'NOT' and instrument == 'NOTCAM' and hdr.get('OBS_MODE', '').lower() == 'imaging':
				instrument_name = 'NOTCAM'
				break

			elif hdr.get('FPA.TELESCOPE') == 'PS1' and hdr.get('FPA.INSTRUMENT') == 'GPC1':
				instrument_name = 'PS1'
				break

			elif telescope == 'Liverpool Telescope':
				instrument_name = 'Liverpool'
				break

			elif telescope == 'CA 3.5m' and instrument == 'Omega2000':
				instrument_name = 'Omega2000'
				break

			elif telescope == 'SWO' and hdr.get('SITENAME') == 'LCO':
				instrument_name = 'Swope'
				break

			elif telescope == 'DUP' and hdr.get('SITENAME') == 'LCO' and instrument == 'Direct/SITe2K-1':
				instrument_name = 'Dupont'
				break

			elif telescope == 'DUP' and instrument == 'RetroCam':
				instrument_name = 'RetroCam'
				break

			elif telescope == 'Baade' and hdr.get('SITENAME') == 'LCO' and instrument == 'FourStar':
				instrument_name = 'Baade'
				break

			elif instrument == 'SOFI' and telescope in ('ESO-NTT', 'other') and (
				origin == 'ESO' or origin.startswith('NOAO-IRAF')):
				instrument_name = 'Sofi'
				break

			elif telescope == 'ESO-NTT' and instrument == 'EFOSC' and (origin == 'ESO' or origin.startswith('NOAO-IRAF')):
				instrument_name = 'EFOSC'
				break

			elif telescope == 'SAI-2.5' and instrument == 'ASTRONIRCAM':
				instrument_name = 'AstroNIRCam'
				break

			elif instrument == 'OMEGACAM' and (origin == 'ESO' or origin.startswith('NOAO-IRAF')):
				instrument_name = 'OmegaCam'
				break

			elif instrument == 'ANDICAM-CCD' and hdr.get('OBSERVAT') == 'CTIO':
				instrument_name = 'AndiCam'
				break

			elif telescope == '1.3m PAIRITEL' and instrument == '2MASS Survey cam':
				instrument_name = 'PairTel'
				break

			elif (origin == 'OAdM' or origin.startswith('NOAO-IRAF')) and telescope == 'TJO' and instrument in (
				'MEIA3', 'MEIA2'):
				instrument_name = 'TJO'
				break

			else:
				raise RuntimeError("Could not determine origin of image")

		image = FlowsImage(image=np.asarray(hdul[ext].data, dtype='float64'), header=hdr, mask=mask)
		ins = instruments[instrument_name](image)

		return ins.image


# get image and wcs solution
# with fits.open(FILENAME, mode='readonly') as hdul:

# --------------------------------------------------------------------------------------------------
def load_image(FILENAME, target_coord=None):
	"""
	Load FITS image.

	Parameters:
		FILENAME (str): Path to FITS file to be loaded.
		target_coord (:class:`astropy.coordinates.SkyCoord`): Coordinates of target.
			Only used for HAWKI images to determine which image extension to load,
			for all other images it is ignored.

	Returns:
		object: Image constainer.

	.. codeauthor:: Rasmus Handberg <rasmush@phys.au.dk>
	"""

	logger = logging.getLogger(__name__)

	# Get image and WCS, find stars, remove galaxies
	image = type('image', (object,), dict())  # image container

	# get image and wcs solution
	with fits.open(FILENAME, mode='readonly') as hdul:

		hdr = hdul[0].header
		image.header = hdr
		origin = hdr.get('ORIGIN', '')
		telescope = hdr.get('TELESCOP', '')
		instrument = hdr.get('INSTRUME', '')

		# Load image data:
		image.image = np.asarray(hdul[0].data, dtype='float64')
		image.shape = image.image.shape

		# Load image mask:
		if origin == 'LCOGT':
			if 'BPM' in hdul:
				image.mask = np.asarray(hdul['BPM'].data, dtype='bool')
			else:
				logger.warning('LCOGT image does not contain bad pixel map. Not applying mask.')
				image.mask = np.zeros_like(image.image, dtype='bool')
		else:
			image.mask = np.zeros_like(image.image, dtype='bool')

		image.mask |= ~np.isfinite(image.image)

		# World Coordinate System:
		with warnings.catch_warnings():
			warnings.simplefilter('ignore', category=FITSFixedWarning)
			image.wcs = WCS(header=hdr, relax=True)

		# Values which will be filled out below, depending on the instrument:
		image.exptime = hdr.get('EXPTIME', None)  # Exposure time * u.second
		image.peakmax = None  # Maximum value above which data is not to be trusted

		# Timestamp:
		if origin == 'LCOGT':
			sites = api.sites.get_all_sites()
			site_keywords = {s['site_keyword']: s for s in sites}
			image.site = site_keywords.get(hdr['SITE'], None)

			observatory = coords.EarthLocation.from_geodetic(lat=hdr['LATITUDE'], lon=hdr['LONGITUD'],
															 height=hdr['HEIGHT'])
			image.obstime = Time(hdr['DATE-OBS'], format='isot', scale='utc', location=observatory)
			image.obstime += 0.5 * image.exptime * u.second  # Make time centre of exposure

			image.photfilter = {
				'zs': 'zp'
			}.get(hdr['FILTER'], hdr['FILTER'])

			# Get non-linear limit
			# TODO: Use actual or some fraction of the non-linearity limit
			# image.peakmax = hdr.get('MAXLIN') # Presumed non-linearity limit from header
			image.peakmax = 60000  # From experience, this one is better.

		elif origin == 'ESO-PARANAL' and telescope == 'ESO-VLT-U4' and instrument == 'HAWKI' and hdr.get(
			'PRODCATG') == 'SCIENCE.MEFIMAGE':
			image.site = api.get_site(2)  # Hard-coded the siteid for ESO Paranal, VLT, UT4
			image.obstime = Time(hdr['DATE-OBS'], format='isot', scale='utc', location=image.site['EarthLocation'])
			image.obstime += 0.5 * image.exptime * u.second  # Make time centre of exposure
			image.photfilter = hdr['FILTER']

			# For HAWKI multi-extension images we search the extensions for which one contains
			# the target, and overwrites the image data with that:
			if target_coord is None:
				raise ValueError("TARGET_COORD is needed for HAWKI images to find the correct extension")
			target_radec = [[target_coord.icrs.ra.deg, target_coord.icrs.dec.deg]]
			for k in range(1, 5):
				w = WCS(header=hdul[k].header, relax=True)
				s = [hdul[k].header['NAXIS2'], hdul[k].header['NAXIS1']]
				pix = w.all_world2pix(target_radec, 0).flatten()
				if pix[0] >= -0.5 and pix[1] >= -0.5 and pix[0] <= s[1] - 0.5 and pix[1] <= s[0] - 0.5:
					ob_type = hdul[k].header["HIERARCH ESO OCS DET1 IMGNAME"].split('_')[-1]
					if "Auto" in ob_type:
						image.ob_type = 'Autojitter'
					elif "Fixed" in ob_type:
						image.ob_type = 'FixedOffset'
					# Should we use a default instead of raising?
					else:
						raise RuntimeError("Image OB Type not AutoJitter or FixedOffset")
					image.image = np.asarray(hdul[k].data, dtype='float64')
					image.shape = image.image.shape
					image.wcs = w
					image.mask = ~np.isfinite(image.image)
					break
			else:
				raise RuntimeError("Could not find image extension that target is on")

		elif telescope == 'NOT' and instrument in ('ALFOSC FASU', 'ALFOSC_FASU') and hdr.get('OBS_MODE',
																							 '').lower() == 'imaging':
			image.site = api.get_site(5)  # Hard-coded the siteid for NOT
			image.obstime = Time(hdr['DATE-AVG'], format='isot', scale='utc', location=image.site['EarthLocation'])

			# Sometimes data from NOT does not have the FILTER keyword,
			# in which case we have to try to figure out which filter
			# was used based on some of the other headers:
			if 'FILTER' in hdr:
				image.photfilter = {
					'B Bes': 'B',
					'V Bes': 'V',
					'R Bes': 'R',
					'g SDSS': 'gp',
					'r SDSS': 'rp',
					'i SDSS': 'ip',
					'i int': 'ip',  # Interference filter
					'u SDSS': 'up',
					'z SDSS': 'zp'
				}.get(hdr['FILTER'].replace('_', ' '), hdr['FILTER'])
			else:
				filters_used = []
				for check_headers in ('ALFLTNM', 'FAFLTNM', 'FBFLTNM'):
					if hdr.get(check_headers) and hdr.get(check_headers).strip().lower() != 'open':
						filters_used.append(hdr.get(check_headers).strip())
				if len(filters_used) == 1:
					image.photfilter = {
						'B_Bes 440_100': 'B',
						'V_Bes 530_80': 'V',
						'R_Bes 650_130': 'R',
						"g'_SDSS 480_145": 'gp',
						"r'_SDSS 618_148": 'rp',
						"i'_SDSS 771_171": 'ip',
						'i_int 797_157': 'ip',  # Interference filter
						"z'_SDSS 832_LP": 'zp'
					}.get(filters_used[0].replace('  ', ' '), filters_used[0])
				else:
					raise RuntimeError("Could not determine filter used.")

			# Get non-linear limit
			# Obtained from http://www.not.iac.es/instruments/detectors/CCD14/LED-linearity/20181026-200-1x1.pdf
			# TODO: grab these from a table for all detector setups of ALFOSC
			image.peakmax = 80000  # For ALFOSC D, 1x1, 200; the standard for SNe.

		elif telescope == 'NOT' and instrument == 'NOTCAM' and hdr.get('OBS_MODE', '').lower() == 'imaging':
			image.site = api.get_site(5)  # Hard-coded the siteid for NOT
			image.obstime = Time(hdr['DATE-AVG'], format='isot', scale='utc', location=image.site['EarthLocation'])

			# Does NOTCAM data sometimes contain a FILTER header?
			# if not we have to try to figure out which filter
			# was used based on some of the other headers:
			if 'FILTER' in hdr:
				raise RuntimeError("NOTCAM: Filter keyword defined")
			filters_used = []
			for check_headers in ('NCFLTNM1', 'NCFLTNM2'):
				if hdr.get(check_headers) and hdr.get(check_headers).strip().lower() != 'open':
					filters_used.append(hdr.get(check_headers).strip())
			if len(filters_used) == 1:
				image.photfilter = {
					'Ks': 'K'
				}.get(filters_used[0], filters_used[0])
			else:
				raise RuntimeError("Could not determine filter used.")

			# Mask out "halo" of pixels with zero value along edge of image:
			image.mask |= edge_mask(image.image, value=0)

		elif hdr.get('FPA.TELESCOPE') == 'PS1' and hdr.get('FPA.INSTRUMENT') == 'GPC1':
			image.site = api.get_site(6)  # Hard-coded the siteid for Pan-STARRS1
			image.obstime = Time(hdr['MJD-OBS'], format='mjd', scale='utc', location=image.site['EarthLocation'])

			image.photfilter = {
				'g.00000': 'gp',
				'r.00000': 'rp',
				'i.00000': 'ip',
				'z.00000': 'zp'
			}.get(hdr['FPA.FILTER'], hdr['FPA.FILTER'])

		elif telescope == 'Liverpool Telescope':
			# Liverpool telescope
			image.site = api.get_site(8)  # Hard-coded the siteid for Liverpool Telescope
			image.obstime = Time(hdr['DATE-OBS'], format='isot', scale='utc', location=image.site['EarthLocation'])
			image.obstime += 0.5 * image.exptime * u.second  # Make time centre of exposure
			image.photfilter = {
				'Bessel-B': 'B',
				'Bessell-B': 'B',
				'Bessel-V': 'V',
				'Bessell-V': 'V',
				'SDSS-U': 'up',
				'SDSS-G': 'gp',
				'SDSS-R': 'rp',
				'SDSS-I': 'ip',
				'SDSS-Z': 'zp'
			}.get(hdr['FILTER1'], hdr['FILTER1'])

		elif telescope == 'CA 3.5m' and instrument == 'Omega2000':
			# Calar Alto 3.5m (Omege2000)
			image.site = api.get_site(9)  # Hard-coded the siteid for Calar Alto 3.5m
			image.obstime = Time(hdr['MJD-OBS'], format='mjd', scale='utc', location=image.site['EarthLocation'])
			image.obstime += 0.5 * image.exptime * u.second  # Make time centre of exposure
			image.photfilter = hdr['FILTER']

		elif telescope == 'SWO' and hdr.get('SITENAME') == 'LCO':
			image.site = api.get_site(10)  # Hard-coded the siteid for Swope, Las Campanas Observatory
			image.obstime = Time(hdr['JD'], format='jd', scale='utc', location=image.site['EarthLocation'])
			image.photfilter = {
				'u': 'up',
				'g': 'gp',
				'r': 'rp',
				'i': 'ip',
			}.get(hdr['FILTER'], hdr['FILTER'])

		elif telescope == 'DUP' and hdr.get('SITENAME') == 'LCO' and instrument == 'Direct/SITe2K-1':
			image.site = api.get_site(14)  # Hard-coded the siteid for Du Pont, Las Campanas Observatory
			image.obstime = Time(hdr['JD'], format='jd', scale='utc', location=image.site['EarthLocation'])
			image.photfilter = {
				'u': 'up',
				'g': 'gp',
				'r': 'rp',
				'i': 'ip',
			}.get(hdr['FILTER'], hdr['FILTER'])

		elif telescope == 'DUP' and instrument == 'RetroCam':
			image.site = api.get_site(16)  # Hard-coded the siteid for Du Pont, Las Campanas Observatory
			image.obstime = Time(hdr['JD'], format='jd', scale='utc', location=image.site['EarthLocation'])
			image.photfilter = {
				'Yc': 'Y',
				'Hc': 'H',
				'Jo': 'J',
			}.get(hdr['FILTER'], hdr['FILTER'])

		elif telescope == 'Baade' and hdr.get('SITENAME') == 'LCO' and instrument == 'FourStar':
			image.site = api.get_site(11)  # Hard-coded the siteid for Swope, Las Campanas Observatory
			image.obstime = Time(hdr['JD'], format='jd', scale='utc', location=image.site['EarthLocation'])
			image.photfilter = {
				'Ks': 'K',
				'J1': 'Y',
			}.get(hdr['FILTER'], hdr['FILTER'])
			image.exptime *= int(hdr['NCOMBINE'])  # EXPTIME is only for a single exposure

		elif instrument == 'SOFI' and telescope in ('ESO-NTT', 'other') and (
			origin == 'ESO' or origin.startswith('NOAO-IRAF')):
			image.site = api.get_site(12)  # Hard-coded the siteid for SOFT, ESO NTT
			if 'TMID' in hdr:
				image.obstime = Time(hdr['TMID'], format='mjd', scale='utc', location=image.site['EarthLocation'])
			else:
				image.obstime = Time(hdr['MJD-OBS'], format='mjd', scale='utc', location=image.site['EarthLocation'])
				image.obstime += 0.5 * image.exptime * u.second  # Make time centre of exposure

			# Photometric filter:
			photfilter_translate = {
				'Ks': 'K'
			}
			if 'FILTER' in hdr:
				image.photfilter = photfilter_translate.get(hdr['FILTER'], hdr['FILTER'])
			else:
				filters_used = []
				for check_headers in ('ESO INS FILT1 ID', 'ESO INS FILT2 ID'):
					if hdr.get(check_headers) and hdr.get(check_headers).strip().lower() != 'open':
						filters_used.append(hdr.get(check_headers).strip())
				if len(filters_used) == 1:
					image.photfilter = photfilter_translate.get(filters_used[0], filters_used[0])
				else:
					raise RuntimeError("Could not determine filter used.")

			# Mask out "halo" of pixels with zero value along edge of image:
			image.mask |= edge_mask(image.image, value=0)

		elif telescope == 'ESO-NTT' and instrument == 'EFOSC' and (origin == 'ESO' or origin.startswith('NOAO-IRAF')):
			image.site = api.get_site(15)  # Hard-coded the siteid for EFOSC, ESO NTT
			image.obstime = Time(hdr['DATE-OBS'], format='isot', scale='utc', location=image.site['EarthLocation'])
			image.obstime += 0.5 * image.exptime * u.second  # Make time centre of exposure
			image.photfilter = {
				'g782': 'gp',
				'r784': 'rp',
				'i705': 'ip',
				'B639': 'B',
				'V641': 'V'
			}.get(hdr['FILTER'], hdr['FILTER'])

		elif telescope == 'SAI-2.5' and instrument == 'ASTRONIRCAM':
			image.site = api.get_site(13)  # Hard-coded the siteid for Caucasus Mountain Observatory
			if 'MIDPOINT' in hdr:
				image.obstime = Time(hdr['MIDPOINT'], format='isot', scale='utc', location=image.site['EarthLocation'])
			else:
				image.obstime = Time(hdr['MJD-AVG'], format='mjd', scale='utc', location=image.site['EarthLocation'])
			image.photfilter = {
				'H_Open': 'H',
				'K_Open': 'K',
			}.get(hdr['FILTER'], hdr['FILTER'])
			image.exptime = hdr.get('FULL_EXP', image.exptime)

		elif instrument == 'OMEGACAM' and (origin == 'ESO' or origin.startswith('NOAO-IRAF')):
			image.site = api.get_site(18)  # Hard-coded the siteid for ESO VLT Survey telescope
			image.obstime = Time(hdr['MJD-OBS'], format='mjd', scale='utc', location=image.site['EarthLocation'])
			image.obstime += 0.5 * image.exptime * u.second  # Make time centre of exposure
			image.photfilter = {
				'i_SDSS': 'ip'
			}.get(hdr['ESO INS FILT1 NAME'], hdr['ESO INS FILT1 NAME'])

		elif instrument == 'ANDICAM-CCD' and hdr.get('OBSERVAT') == 'CTIO':
			image.site = api.get_site(
				20)  # Hard-coded the siteid for ANDICAM at Cerro Tololo Interamerican Observatory (CTIO)
			image.obstime = Time(hdr['JD'], format='jd', scale='utc', location=image.site['EarthLocation'])
			image.obstime += 0.5 * image.exptime * u.second  # Make time centre of exposure
			image.photfilter = hdr['CCDFLTID']

		elif telescope == '1.3m PAIRITEL' and instrument == '2MASS Survey cam':
			image.site = api.get_site(21)  # Hard-coded the siteid for Peters Automated InfraRed Imaging TELescope
			time_start = Time(hdr['STRT_CPU'], format='iso', scale='utc', location=image.site['EarthLocation'])
			time_stop = Time(hdr['STOP_CPU'], format='iso', scale='utc', location=image.site['EarthLocation'])
			image.obstime = time_start + 0.5 * (time_stop - time_start)
			image.photfilter = {
				'j': 'J',
				'h': 'H',
				'k': 'K',
			}.get(hdr['FILTER'], hdr['FILTER'])

			# Mask out "halo" of pixels with zero value along edge of image:
			image.mask |= edge_mask(image.image, value=0)

		elif (origin == 'OAdM' or origin.startswith('NOAO-IRAF')) and telescope == 'TJO' and instrument in (
			'MEIA3', 'MEIA2'):
			image.site = api.get_site(
				22)  # Hard-coded the siteid for Telescopi Joan Oró (TJO) at Observatori Astronòmic del Montsec
			image.obstime = Time(hdr['JD'], format='jd', scale='utc', location=image.site['EarthLocation'])
			image.obstime += 0.5 * image.exptime * u.second  # Make time centre of exposure
			image.photfilter = hdr['FILTER']

		else:
			raise RuntimeError("Could not determine origin of image")

	# Sanity checks:
	if image.exptime is None:
		raise ValueError("Image exposure time could not be extracted")

	# Create masked version of image:
	image.image[image.mask] = np.NaN
	image.clean = np.ma.masked_array(data=image.image, mask=image.mask, copy=False)

	return image


def main():
	filename = "../ADP.2021-11-08T12_34_46.383.fits"
	img = new_load_image(filename, (325.772590825, -17.54410925))
	print(img.peakmax)
	return img


if __name__ == '__main__':
	main()
