from typing import Optional, Callable

import numpy as np
from astropy.stats import SigmaClip
from numpy.typing import ArrayLike
from photutils.background import Background2D, SExtractorBackground
from photutils.utils import calc_total_error


class FlowsBackground:

    def __init__(self, background_estimator: Background2D = Background2D):
        self.background_estimator = Background2D
        self.background: Optional[ArrayLike] = None
        self.background_rms: Optional[ArrayLike] = None

    def estimate_background(self, clean_image: np.ma.MaskedArray) -> None:
        # Estimate image background:
        # Not using image.clean here, since we are redefining the mask anyway
        bkg2d = self.background_estimator(clean_image, (128, 128), filter_size=(5, 5),
                                         sigma_clip=SigmaClip(sigma=3.0), bkg_estimator=SExtractorBackground(),
                                         exclude_percentile=50.0)
        self.background = bkg2d.background
        self.background_rms = bkg2d.background_rms

    def background_subtract(self, clean_image: ArrayLike) -> ArrayLike:
        if self.background is None:
            self.estimate_background(clean_image)
        return clean_image - self.background

    def error(self, clean_image: ArrayLike, error_method: Callable = calc_total_error):
        """
        Calculate the 2D error using the background RMS.
        """
        if self.background is None:
            raise AttributeError("background must be estimated before calling error")
        return error_method(clean_image, self.background_rms, 1.0)
