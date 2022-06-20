"""
Load image code.
"""
from __future__ import annotations
import numpy as np
import logging
from typing import Union, Tuple
import astropy.coordinates as coords
from astropy.io import fits
from astropy.time import Time

from .instruments import INSTRUMENTS
from .image import FlowsImage

logger = logging.getLogger(__name__)  # Singleton logger instance


def load_image(filename: str, target_coord: Union[coords.SkyCoord, Tuple[float, float]] = None):
    """
    Load FITS image using FlowsImage class and Instrument Classes.

    Parameters:
        filename (str): Path to FITS file to be loaded.
        target_coord (:class:`astropy.coordinates.SkyCoord`): Coordinates of target.
            Only used for HAWKI images to determine which image extension to load,
            for all other images it is ignored.

    Returns:
        FlowsImage: instance of FlowsImage with values populated based on instrument.

    """
    ext = 0  # Default extension is  0, individual instruments may override this.
    # Read fits image, Structural Pattern Match to specific instrument.
    with fits.open(filename, mode='readonly') as hdul:
        hdr = hdul[ext].header
        origin = hdr.get('ORIGIN', '')
        telescope = hdr.get('TELESCOP', '')
        instrument = hdr.get('INSTRUME', '')

        for inst_name, inst_cls in INSTRUMENTS:
            if inst_cls.identifier(telescope, origin, instrument, hdr):
                logger.info(f"Image is using instrument {inst_name}")
                ext = inst_cls.get_ext(hdul, target_coord)
                mask = inst_cls.get_mask(hdul)
                # Default = None is to only mask all non-finite values, override here is additive.

                image = FlowsImage(image=np.asarray(hdul[ext].data, dtype='float64'),
                                   header=hdr, mask=mask)
                current_instrument = inst_cls(image)
                clean_image = current_instrument.process_image()
                if target_coord is not None:
                    clean_image.obstime = correct_barycentric(clean_image.obstime, target_coord)
                return clean_image

        raise RuntimeError("Could not determine origin of image")

def correct_barycentric(obstime: Time, target_coord: coords.SkyCoord) -> Time:
    """
    BARYCENTRIC CORRECTION OF TIME

    Parameters:
        obstime (astropy.time.Time): Midpoint observed image time.
        target_coord (astropy.coords.SkyCoord): Coordinates of target in image.

    Returns:
        obstime (astropy.time.Time): Time corrected to barycenter with jpl ephemeris.
    """
    ltt_bary = obstime.light_travel_time(target_coord, ephemeris='jpl')
    return obstime.tdb + ltt_bary
