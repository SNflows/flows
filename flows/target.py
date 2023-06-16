from dataclasses import dataclass
from typing import Dict, Optional

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from numpy.typing import NDArray
from tendrils import api


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
                               pixel_pos: Optional[NDArray] = None) -> None:
        """
        Add pixel coordinates to target.
        """
        if pixel_column is None or pixel_row is None:
            if pixel_pos is None:
                raise ValueError('Either pixel_column, pixel_row or pixel_pos must be provided.')
            pixel_column, pixel_row = pixel_pos

        self.pixel_column = pixel_column
        self.pixel_row = pixel_row

    def output_dict(self, starid: Optional[int] = 0) -> Dict:
        """
        Return target as output dictionary. starid = -1 means difference image.
        """
        return {'starid': starid, 'ra': self.ra, 'decl': self.dec, 'pixel_column': self.pixel_column,
                'pixel_row': self.pixel_row, 'used_for_epsf': False}

    @classmethod
    def from_dict(cls, d: Dict) -> 'Target':
        """
        Create target from dictionary.
        """
        return cls(ra=d['ra'], dec=d['decl'], name=d['target_name'], id=d['targetid'], photfilter=d['photfilter'],)

    @classmethod
    def from_fid(cls, fid: int, datafile: Optional[Dict] = None) -> 'Target':
        """
        Create target from fileid.
        """

        datafile = datafile or api.get_datafile(fid)
        if datafile is None:
            raise ValueError(f'No datafile found for fid={fid}')
        d = api.get_target(datafile['target_name']) | datafile
        return cls.from_dict(d)

    @classmethod
    def from_tid(cls, target_id: int) -> 'Target':
        """
        Create target from target id.
        """
        target_pars = api.get_target(target_id)
        return cls(
            ra=target_pars['ra'], dec=target_pars['decl'],
            name=target_pars['target_name'], id=target_pars['targetid'])
