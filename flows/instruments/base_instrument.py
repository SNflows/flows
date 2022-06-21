"""
Base instrument class
"""
# Standard lib
from __future__ import annotations
import logging
from typing import Union, Dict, Any, Optional
from abc import ABC, abstractmethod
# Third party
import astropy.coordinates as coords
from astropy.io import fits
from astropy.time import Time
# First party
from tendrils import api
from flows.image import FlowsImage, InstrumentDefaults
from flows.utilities import create_logger
logger = create_logger()


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

    @abstractmethod
    def process_image(self):
        pass


class Instrument(AbstractInstrument):
    peakmax: int = None
    siteid: int = None
    telescope: str = ''  # Fits Header name of TELESCOP
    instrument: str = ''  # Fits Header name of Instrument (can be partial)
    origin: str = ''  # Fits Header value of ORIGIN (if relevant)
    unique_headers: Optional[Dict[str, Any]] = None  # Unique key value pairs from header for identifying instrument.

    def __init__(self, image: FlowsImage = None, header: fits.header.Header = None):
        self.image = image
        self.hdr = header

    def get_site(self) -> Dict[str, Any]:
        if self.siteid is not None:
            return api.get_site(self.siteid)

    def get_exptime(self) -> Union[float, int, str]:
        exptime = self.image.header.get('EXPTIME', None)
        if exptime is None:
            raise ValueError("Image exposure time could not be extracted")
        return exptime

    def get_obstime(self) -> Time:
        """Default for JD, jd, utc."""
        return Time(self.image.header['JD'], format='jd', scale='utc', location=self.image.site['EarthLocation'])

    def get_photfilter(self):
        return self.image.header['FILTER']

    def set_instrument_defaults(self):
        """
        Set default values for instrument.
        """
        self.image.instrument_defaults = InstrumentDefaults()

    def _get_clean_image(self):
        self.image.peakmax = self.peakmax
        self.image.site = self.get_site()
        self.image.exptime = self.get_exptime()
        self.image.obstime = self.get_obstime()
        self.image.photfilter = self.get_photfilter()
        self.image.create_masked_image()

    def process_image(self, image: FlowsImage = None) -> FlowsImage:
        """Process existing or new image."""
        if image is not None:
            self.image = image
        if self.image is None:
            raise AttributeError('No FlowsImage to be processed. Self.image was None')

        self._get_clean_image()
        self.set_instrument_defaults()
        return self.image

    @classmethod
    def identifier(cls, telescope: str, origin: str, instrument: str, hdr: fits.header.Header) -> bool:
        """Unique identifier"""
        unique_conds = all([hdr.get(key) == cls.unique_headers.get(key) for key in cls.unique_headers.keys()] if
                           cls.unique_headers is not None else [True])

        return all([cls.telescope in telescope if cls.telescope != '' else True,
                    cls.origin == origin if cls.origin != '' else True,
                    cls.instrument in instrument if cls.instrument != '' else True,
                    unique_conds])

    @staticmethod
    def get_ext(hdul: fits.HDUList, target_coords: coords.SkyCoord = None) -> int:
        """Instruments which need a special treatment to find the image extension
        should overwrite this."""
        return 0

    @staticmethod
    def get_mask(hdul: fits.HDUList) -> Optional[Any]:
        """Instruments which need a special treatment to find the mask should overwrite this."""
        return None
