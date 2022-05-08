from dataclasses import dataclass
from typing import Optional, Dict

import numpy as np
from numpy.typing import ArrayLike, NDArray
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS


@dataclass
class Target:
    ra: float
    dec: float
    name: Optional[str] = None
    id: Optional[int] = None  # Target id from Flows database
    photfilter: Optional[str] = None  # Defined if target is associated with an image.
    coords: Optional[SkyCoord] = None
    pixel_column: Optional[int] = None
    pixel_row: Optional[int] = None

    def __post_init__(self):
        if self.coords is None:
            self.coords = SkyCoord(ra=self.ra, dec=self.dec, unit='deg', frame='icrs')

    def calc_pixels(self, wcs: WCS) -> None:
        pixels = np.array(wcs.all_world2pix(self.ra, self.dec, 1)).T
        self._add_pixel_coordinates(pixel_pos=pixels)

    def _add_pixel_coordinates(self, pixel_column: Optional[int] = None, pixel_row: Optional[int] = None,
                               pixel_pos: Optional[ArrayLike] = None) -> None:
        """
        Add pixel coordinates to target.
        """
        if pixel_pos:
            self.pixel_column = pixel_pos[0]
            self.pixel_row = pixel_pos[1]
        elif pixel_column is not None and pixel_row is not None:
            self.pixel_column = pixel_column
            self.pixel_row = pixel_row
        raise AttributeError("Either pixel_pos or pixel_column and pixel_row must be defined.")

    def output_dict(self, starid: Optional[int] = 0) -> Dict:
        """
        Return target as output dictionary. starid = -1 means difference image.
        """
        return {'starid': starid, 'ra': self.ra, 'decl': self.dec, 'pixel_column': self.pixel_column,
                'pixel_row': self.pixel_row, 'used_for_epsf': False}
