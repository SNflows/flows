from __future__ import annotations
import numpy as np
import logging
from dataclasses import dataclass  # , field
import warnings
from typing import Union
from astropy.time import Time
from astropy.wcs import WCS, FITSFixedWarning
from typing import Tuple,  Dict, Any, Optional


logger = logging.getLogger(__name__)  # Singleton logger instance


@dataclass
class InstrumentDefaults:
    """
    Default radius and FWHM for an instrument in arcseconds.
    """
    radius: float = 10
    fwhm: float = 6.0   # Best initial guess
    fwhm_min: float = 3.5
    fwhm_max: float = 18.0


@dataclass
class FlowsImage:
    image: np.ndarray
    header: Dict
    mask: Optional[np.ndarray] = None
    peakmax: Optional[float] = None
    exptime: Optional[float] = None
    instrument_defaults: Optional[InstrumentDefaults] = None
    site: Optional[Dict[str, Any]] = None
    obstime: Optional[Time] = None
    photfilter: Optional[str] = None
    wcs: Optional[WCS] = None

    clean: Optional[np.ma.MaskedArray] = None
    subclean: Optional[np.ma.MaskedArray] = None
    error: Optional[np.ma.MaskedArray] = None

    def __post_init__(self):
        self.shape = self.image.shape
        self.wcs = self.create_wcs()
        # Make empty mask
        if self.mask is None:
            self.mask = np.zeros_like(self.image, dtype='bool')
        self.check_finite()

    def check_finite(self):
        self.mask |= ~np.isfinite(self.image)

    def update_mask(self, mask):
        self.mask = mask
        self.check_finite()

    def create_wcs(self) -> WCS:
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', category=FITSFixedWarning)
            return WCS(header=self.header, relax=True)

    def create_masked_image(self):
        """Warning: this is destructive and will overwrite image data setting masked values to NaN"""
        self.image[self.mask] = np.NaN
        self.clean = np.ma.masked_array(data=self.image, mask=self.mask, copy=False)

    def set_edge_rows_to_value(self, y: Tuple[float] = None, value: Union[int, float, np.float64] = 0):
        if y is None:
            pass
        for row in y:
            self.image[row] = value

    def set_edge_columns_to_value(self, x: Tuple[float] = None, value: Union[int, float, np.float64] = 0):
        if x is None:
            pass
        for col in x:
            self.image[:, col] = value

    @staticmethod
    def get_edge_mask(img: np.ndarray, value: Union[int, float, np.float64] = 0):
        """
        Create boolean mask of given value near edge of image.

        Parameters:
            img (ndarray): image with values for masking.
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

    def apply_edge_mask(self, y: Tuple[int] = None, x: Tuple[int] = None, apply_existing_mask_first: bool = False):
        """
        Masks given rows and columns of image but will replace the current mask! Set apply_existing_mask_first to True
        if the current mask should be kept.
        :param y: Tuple[int] of rows to mask
        :param x: Tuple[int] of columns to mask
        :param apply_existing_mask_first: Whether to apply the existing mask to image first, before overwriting mask.
        :return: None
        """
        if y is None and x is None:
            logger.debug("(y,x) was None when applying edge mask. Edge was not actually masked.")

        if apply_existing_mask_first:
            self.create_masked_image()

        if y is not None:
            self.set_edge_rows_to_value(y=y)

        if x is not None:
            self.set_edge_columns_to_value(x=x)

        self.mask = self.get_edge_mask(self.image)
        self.create_masked_image()
