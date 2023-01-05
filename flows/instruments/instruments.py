"""
Instrument classes that inherit from the base class.
Modify to add your own instrument.

Identifying relevant image properties:

`site` = required, add manually.
`peakmax` = optional, add manually or provide header.
The rest of `exptime`, `obstime`, `photfilter` use the
following (overrideable) base functions. Override
if the default one from the baseclass does not fit
your instrument. See:
    ```
    self.image.peakmax = self.peakmax
    self.image.site = self.get_site() 
    self.image.exptime = self.get_exptime()
    self.image.obstime = self.get_obstime()
    self.image.photfilter = self.get_photfilter()
    ```

Identifying the instrument for an image:

Each instrument can define (one or many) of `origin`, 
`telescope`, `instrument` fields correspinding to the 
standard fits headers to help uniquely identify itself.
More advanced logic is possible using `unique_headers`
field as a dict of key,value pairs in the header. Ex:
unique_headers = {'PRODCATG': 'SCIENCE.MEFIMAGE'}.
These are all optional, defaults are set in baseclass.
"""
# Standard lib
from __future__ import annotations

import inspect
import sys
from typing import List, Optional, Tuple, Union

import astropy.coordinates as coords
import astropy.units as u
# Third party
import numpy as np
from astropy.io import fits
from astropy.time import Time
from astropy.wcs import WCS
# First party
from tendrils import api

from flows.filters import FILTERS
from flows.image import FlowsImage
from flows.instruments.base_instrument import Instrument
from flows.utilities import create_logger

logger = create_logger()


class LCOGT(Instrument):
    siteid = 1  # Can be between 1, 3, 4, 6, 17, 19. @TODO: Refactor to own classes.
    peakmax: int = 60000
    origin = 'LCOGT'

    def get_site(self):
        nonesite = {'siteid': None}
        if self.image is None:
            return nonesite
        sites = api.sites.get_all_sites()
        site_keywords = {s['site_keyword']: s for s in sites}
        site = site_keywords.get(self.image.header['SITE'], nonesite)
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

    @staticmethod
    def get_mask(hdul: fits.HDUList):
        if 'BPM' in hdul:
            return np.asarray(hdul['BPM'].data, dtype='bool')

        logger.warning('LCOGT image does not contain bad pixel map. Not applying mask.')
        return None


def verify_coordinates(target_coords: Union[coords.SkyCoord, Tuple]) -> Optional[coords.SkyCoord]:
    if target_coords is None:
        return None
    if isinstance(target_coords, coords.SkyCoord):
        return target_coords
    elif len(target_coords) == 2:
        return coords.SkyCoord(ra=target_coords[0] * u.deg, dec=target_coords[1] * u.deg, frame='icrs')
    return None


class HAWKI(Instrument):
    siteid = 2  # Hard-coded the siteid for ESO Paranal, VLT, UT4
    telescope = 'ESO-VLT-U4'  # Fits Header name of TELESCOP
    instrument = 'HAWKI'  # Fits Header name of Instrument (can be partial)
    origin = 'ESO-PARANAL'  # Fits Header value of ORIGIN (if relevant)
    #unique_headers = {'PRODCATG': 'SCIENCE.MEFIMAGE'}

    def __init__(self, image: FlowsImage = None):
        super().__init__(image)
        if self.image is not None:
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

    @staticmethod
    def get_ext(hdul: fits.HDUList, target_coords: coords.SkyCoord = None,
                fallback_extension: int = None) -> int:
        target_coord = verify_coordinates(target_coords)
        if target_coord is None:
            raise ValueError("TARGET_COORD is needed for HAWKI images to find the correct extension")

        # Incase this is not a multi-extension imageL
        if len(hdul) == 1:
            return 0

        # For HAWKI multi-extension images we search the extensions for which one contains
        # the target, Create Image from that extension.
        target_radec = [[target_coord.icrs.ra.deg, target_coord.icrs.dec.deg]]

        for k in range(1, 5):
            w = WCS(header=hdul[k].header, relax=True)
            s = [hdul[k].header['NAXIS2'], hdul[k].header['NAXIS1']]
            pix = w.all_world2pix(target_radec, 0).flatten()
            if -0.5 <= pix[0] <= s[1] - 0.5 and -0.5 <= pix[1] <= s[0] - 0.5:
                return k
        if fallback_extension is not None:
            return fallback_extension
        else:
            raise RuntimeError(f"Could not find image extension that target is on!")



class ALFOSC(Instrument):
    # Obtained from http://www.not.iac.es/instruments/detectors/CCD14/LED-linearity/20181026-200-1x1.pdf
    peakmax = 80000  # For ALFOSC D, 1x1, 200; the standard for SNe.
    siteid = 5
    telescope = "NOT"
    instrument = "ALFOSC"
    unique_headers = {"OBS_MODE": 'IMAGING'}

    def get_obstime(self):
        return Time(self.image.header['DATE-AVG'], format='isot', scale='utc',
                    location=self.image.site['EarthLocation'])

    def get_photfilter(self):
        # Sometimes data from NOT does not have the FILTER keyword,
        # in which case we have to try to figure out which filter
        # was used based on some of the other headers:
        if 'FILTER' in self.image.header:
            photfilter = {'B Bes': 'B', 'V Bes': 'V', 'R Bes': 'R', 'g SDSS': 'gp', 'r SDSS': 'rp', 'i SDSS': 'ip',
                          'i int': 'ip',  # Interference filter
                          'u SDSS': 'up', 'z SDSS': 'zp'}.get(self.image.header['FILTER'].replace('_', ' '),
                                                              self.image.header['FILTER'])
        else:
            filters_used = []
            for check_headers in ('ALFLTNM', 'FAFLTNM', 'FBFLTNM'):
                isopen = self.image.header.get(check_headers).strip().lower() != 'open'
                if self.image.header.get(check_headers) and isopen:
                    filters_used.append(self.image.header.get(check_headers).strip())
            if len(filters_used) == 1:
                photfilter = {'B_Bes 440_100': 'B', 'V_Bes 530_80': 'V', 'R_Bes 650_130': 'R', "g'_SDSS 480_145": 'gp',
                              "r'_SDSS 618_148": 'rp', "i'_SDSS 771_171": 'ip', 'i_int 797_157': 'ip',
                              # Interference filter
                              "z'_SDSS 832_LP": 'zp'}.get(filters_used[0].replace('  ', ' '), filters_used[0])
            else:
                raise RuntimeError("Could not determine filter used.")

        return photfilter


class NOTCAM(Instrument):
    siteid = 5
    telescope = "NOT"
    instrument = "NOTCAM"
    # unique_headers = {"OBS_MODE": 'IMAGING'}  # not needed.

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
            photfilter = {'Ks': 'K'}.get(filters_used[0], filters_used[0])
        else:
            raise RuntimeError("Could not determine filter used.")
        return photfilter


class PS1(Instrument):
    siteid = 6
    unique_headers = {'FPA.TELESCOPE': 'PS1', 'FPA.INSTRUMENT': 'GPC1'}

    def get_obstime(self):
        return Time(self.image.header['MJD-OBS'], format='mjd', scale='utc', location=self.image.site['EarthLocation'])

    def get_photfilter(self):
        photfilter = {'g.00000': 'gp', 'r.00000': 'rp', 'i.00000': 'ip', 'z.00000': 'zp'}.get(
            self.image.header['FPA.FILTER'], self.image.header['FPA.FILTER'])
        return photfilter


class Liverpool(Instrument):
    siteid = 8
    telescope = 'Liverpool Telescope'

    def get_obstime(self):
        obstime = Time(self.image.header['DATE-OBS'], format='isot', scale='utc',
                       location=self.image.site['EarthLocation'])
        obstime += 0.5 * self.image.exptime * u.second  # Make time centre of exposure
        return obstime

    def get_photfilter(self):
        photfilter = {'Bessel-B': 'B', 'Bessell-B': 'B', 'Bessel-V': 'V', 'Bessell-V': 'V', 'SDSS-U': 'up',
                      'SDSS-G': 'gp', 'SDSS-R': 'rp', 'SDSS-I': 'ip', 'SDSS-Z': 'zp'}.get(self.image.header['FILTER1'],
                                                                                          self.image.header['FILTER1'])
        return photfilter


class Omega2000(Instrument):
    siteid = 9
    telescope = 'CA 3.5m'
    instrument = 'Omega2000'

    def get_obstime(self):
        obstime = Time(self.image.header['MJD-OBS'], format='mjd', scale='utc',
                       location=self.image.site['EarthLocation'])
        obstime += 0.5 * self.image.exptime * u.second
        return obstime


class Swope(Instrument):
    siteid = 10
    telescope = "SWO"

    def get_photfilter(self):
        photfilter = {'u': 'up', 'g': 'gp', 'r': 'rp', 'i': 'ip', }.get(self.image.header['FILTER'],
                                                                        self.image.header['FILTER'])
        return photfilter

    @classmethod
    def identifier(cls, telescope: str, origin: str, instrument: str, hdr: fits.header.Header) -> bool:
        """Unique identifier"""
        return telescope.upper().startswith(cls.telescope) and hdr.get('SITENAME') == 'LCO'


class Swope_newheader(Swope):

    def get_obstime(self):
        obstime = Time(self.image.header['MJD-OBS'], format='mjd', scale='utc',
                       location=self.image.site['EarthLocation'])
        obstime += 0.5 * self.image.exptime * u.second
        return obstime

    @classmethod
    def identifier(cls, telescope: str, origin: str, instrument: str, hdr: fits.header.Header) -> bool:
        """Unique identifier"""
        return telescope.upper().startswith('SWO') and origin == 'ziggy'


class Dupont(Instrument):
    siteid = 14
    telescope = 'DUP'
    instrument = 'Direct/SITe2K-1'
    unique_headers = {'SITENAME': 'LCO'}

    def get_photfilter(self):
        photfilter = {'u': 'up', 'g': 'gp', 'r': 'rp', 'i': 'ip', }.get(self.image.header['FILTER'],
                                                                        self.image.header['FILTER'])
        return photfilter


class RetroCam(Instrument):
    siteid = 16
    telescope = 'DUP'
    instrument = 'RetroCam'

    def get_photfilter(self):
        photfilter = {'Yc': 'Y', 'Hc': 'H', 'Jo': 'J', }.get(self.image.header['FILTER'], self.image.header['FILTER'])
        return photfilter


class Baade(Instrument):
    siteid = 11
    telescope = 'Baade'
    instrument = 'FourStar'
    unique_headers = {'SITENAME': 'LCO'}

    def get_exptime(self):
        exptime = super().get_exptime()
        exptime *= int(self.image.header['NCOMBINE'])  # EXPTIME is only for a single exposure
        return exptime

    def get_photfilter(self):
        photfilter = {'Ks': 'K', 'J1': 'Y', }.get(self.image.header['FILTER'], self.image.header['FILTER'])
        return photfilter


class Sofi(Instrument):
    siteid = 12
    instrument = 'SOFI'

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
        photfilter_translate = {'Ks': 'K'}
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

    @classmethod
    def identifier(cls, telescope, origin, instrument, hdr):
        return instrument == cls.instrument and telescope in ('ESO-NTT', 'other')


class EFOSC(Instrument):
    siteid = 15
    telescope = 'ESO-NTT'
    instrument = 'EFOSC'

    def get_obstime(self):
        obstime = Time(self.image.header['DATE-OBS'], format='isot', scale='utc',
                       location=self.image.site['EarthLocation'])
        obstime += 0.5 * self.image.exptime * u.second  # Make time centre of exposure
        return obstime

    def get_photfilter(self):
        hdr = self.image.header
        photfilter = {'g782': 'gp', 'r784': 'rp', 'i705': 'ip', 'B639': 'B', 'V641': 'V'}.get(hdr['FILTER'],
                                                                                              hdr['FILTER'])
        return photfilter


class AstroNIRCam(Instrument):
    siteid = 13
    telescope = 'SAI-2.5'
    instrument = 'ASTRONIRCAM'

    def get_exptime(self):
        exptime = self.image.header.get('FULL_EXP', None)
        if exptime is not None:
            return exptime
        return super().get_exptime()

    def get_obstime(self):
        hdr = self.image.header
        if 'MIDPOINT' in hdr:
            obstime = Time(hdr['MIDPOINT'], format='isot', scale='utc', location=self.image.site['EarthLocation'])
        else:
            obstime = Time(hdr['MJD-AVG'], format='mjd', scale='utc', location=self.image.site['EarthLocation'])
        return obstime

    def get_photfilter(self):
        hdr = self.image.header
        photfilter = {'H_Open': 'H', 'K_Open': 'K', 'J_Open': 'J'}.get(hdr['FILTER'], hdr['FILTER'])
        return photfilter


class OmegaCam(Instrument):
    siteid = 18  # Hard-coded the siteid for ESO VLT Survey telescope
    instrument = 'OMEGACAM'

    def get_obstime(self):
        obstime = Time(self.image.header['MJD-OBS'], format='mjd', scale='utc',
                       location=self.image.site['EarthLocation'])
        obstime += 0.5 * self.image.exptime * u.second  # Make time centre of exposure
        return obstime

    def get_photfilter(self):
        hdr = self.image.header
        photfilter = {'i_SDSS': 'ip'}.get(hdr['ESO INS FILT1 NAME'], hdr['ESO INS FILT1 NAME'])
        return photfilter


class AndiCam(Instrument):
    siteid = 20  # Hard-coded the siteid for ANDICAM at Cerro Tololo Interamerican Observatory (CTIO)
    instrument = 'ANDICAM-CCD'
    unique_headers = {'OBSERVAT': 'CTIO'}

    def get_obstime(self):
        obstime = super().get_obstime()
        obstime += 0.5 * self.image.exptime * u.second
        return obstime

    def get_photfilter(self):
        return self.image.header['CCDFLTID']


class PairTel(Instrument):
    siteid = 21
    telescope = "1.3m PAIRITEL"
    instrument = "2MASS Survey cam"

    def get_obstime(self):
        hdr = self.image.header
        time_start = Time(hdr['STRT_CPU'], format='iso', scale='utc', location=self.image.site['EarthLocation'])
        time_stop = Time(hdr['STOP_CPU'], format='iso', scale='utc', location=self.image.site['EarthLocation'])
        obstime = time_start + 0.5 * (time_stop - time_start)
        return obstime

    def get_photfilter(self):
        hdr = self.image.header
        photfilter = {'j': 'J', 'h': 'H', 'k': 'K', }.get(hdr['FILTER'], hdr['FILTER'])
        return photfilter


class TJO_MEIA2(Instrument):
    siteid = 22
    telescope = 'TJO'
    instrument = 'MEIA2'


    def get_obstime(self):
        obstime = super().get_obstime()
        obstime += 0.5 * self.image.exptime * u.second
        return obstime


class TJO_MEIA3(Instrument):
    siteid = 22
    telescope = 'TJO'
    instrument = 'MEIA3'

    def get_obstime(self):
        obstime = super().get_obstime()
        obstime += 0.5 * self.image.exptime * u.second
        return obstime


class RATIR(Instrument):
    siteid = 23
    telescope = "OAN/SPM Harold L. Johnson 1.5-meter"

    def get_obstime(self):
        obstime = Time(self.image.header['DATE-OBS'], format='isot', scale='utc',
                        location=self.image.site['EarthLocation'])
        return obstime

    def get_photfilter(self):
        ratir_filt = self.image.header['FILTER']
        if ratir_filt in ['Z', 'r', 'i']:
            return {'Z': 'zp', 'r': 'rp', 'i': 'ip'}.get(ratir_filt)
        return ratir_filt


class AFOSC(Instrument):
    siteid = 25
    peakmax = 50_000
    telescope = '1.82m Reflector'  # Fits Header name of TELESCOP
    instrument = 'AFOSC'  # Fits Header name of Instrument (can be partial)

    def get_obstime(self):
        obstime = Time(self.image.header['DATE-OBS'], format='isot', scale='utc',
                        location=self.image.site['EarthLocation'])
        obstime += 0.5 * self.image.exptime * u.second  # Make time centre of exposure
        return obstime

    def get_photfilter(self):
        filt = self.image.header['FILTER']

        if "sloan" in filt.lower():
            return filt[0]+'p'  # Return gp,rp,ip,zp for g-sloan, etc.
        elif filt in FILTERS.keys():
            return filt
        elif filt+'p' in FILTERS.keys():
            return filt+'p'
        elif filt[0] in FILTERS.keys():
            return filt[0]

        raise ValueError(f"Could not find filter {filt} in {[f for f in FILTERS.keys()]}")


class Schmidt(Instrument):
    siteid = 26
    peakmax = 56_000
    telescope = '67/91 Schmidt Telescope' 
    instrument = 'Moravian G4-16000LC'
    origin = '' 

    def get_obstime(self):
        obstime = Time(self.image.header['DATE-OBS'], format='isot', scale='utc',
                        location=self.image.site['EarthLocation'])
        obstime += 0.5 * self.image.exptime * u.second  # Make time centre of exposure
        return obstime

    def get_photfilter(self):
        filt = self.image.header['FILTER']

        if "sloan" in filt.lower():
            return filt[0]+'p'  # Return gp,rp,ip,zp for g-sloan, etc.
        elif filt in FILTERS.keys():
            return filt
        elif filt+'p' in FILTERS.keys():
            return filt+'p'

        raise ValueError(f"Could not find filter {filt} in {[f for f in FILTERS.keys()]}")


class TNG(Instrument):
    siteid = 5  # same as NOT
    #instrument = 'LRS'
    #telescope = 'TNG'
    
    unique_headers = {'TELESCOP':'TNG', 'INSTRUME':'LRS'} # assume we use unique headers?


    def get_obstime(self):
        return Time(self.image.header['DATE-OBS'], format='isot', scale='utc',
                    location=self.image.site['EarthLocation'])


    def get_exptime(self):
        exptime = super().get_exptime()
        exptime *= int(self.image.header['EXPTIME'])  
        return exptime

    def get_photfilter(self):
        ratir_filt = self.image.header['FLT_ID']
        if ratir_filt in ['B_John_10', 'g_sdss_30', 'r_sdss_31', 'i_sdss_32',
        'u_sdss_29', 'V_John_11']:
            return {'B_John_10': 'B', 'g_sdss_30':'g', 'r':'r_sdss_31',
            'i_sdss_32':'i', 'u_sdss_29':'u', 'V_John_11':'V_John_11' }.get(ratir_filt)
        return ratir_filt


INSTRUMENTS: list[tuple[str, Instrument]] = inspect.getmembers(sys.modules[__name__],
                                 lambda member: inspect.isclass(member) and member.__module__ == __name__)
                                 