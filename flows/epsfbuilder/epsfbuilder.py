#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Photutils hack for EPSF building

.. codeauthor:: Simon Holmbo <simonholmbo@phys.au.dk>
"""
import time
import numpy as np
from scipy.interpolate import griddata, UnivariateSpline
import photutils.psf
from typing import List, Tuple
from utilities import create_logger
logger = create_logger()


class FlowsEPSFBuilder(photutils.psf.EPSFBuilder):
    def _create_initial_epsf(self, stars):
        epsf = super()._create_initial_epsf(stars)
        epsf.origin = None

        X, Y = np.meshgrid(*map(np.arange, epsf.shape[::-1]))

        X = X / epsf.oversampling[0] - epsf.x_origin
        Y = Y / epsf.oversampling[1] - epsf.y_origin

        self._epsf_xy_grid = X, Y

        return epsf

    def _resample_residual(self, star, epsf):
        # max_dist = .5 / np.sqrt(np.sum(np.power(epsf.oversampling, 2)))

        # star_points = list(zip(star._xidx_centered, star._yidx_centered))
        # epsf_points = list(zip(*map(np.ravel, self._epsf_xy_grid)))

        # star_tree = cKDTree(star_points)
        # dd, ii = star_tree.query(epsf_points, distance_upper_bound=max_dist)
        # mask = np.isfinite(dd)

        # star_data = np.full_like(epsf.data, np.nan)
        # star_data.ravel()[mask] = star._data_values_normalized[ii[mask]]

        star_points = list(zip(star._xidx_centered, star._yidx_centered))
        star_data = griddata(star_points, star._data_values_normalized, self._epsf_xy_grid)

        return star_data - epsf._data

    def __call__(self, *args, **kwargs):
        t0 = time.time()

        epsf, stars = super().__call__(*args, **kwargs)

        epsf.fit_info = dict(n_iter=len(self._epsf), max_iters=self.maxiters, time=time.time() - t0, )

        return epsf, stars


def verify_epsf(epsf: photutils.psf.EPSFModel) -> Tuple[bool, List[float]]:
    fwhms = []
    epsf_ok = True
    for a in (0, 1):
        # Collapse the PDF along this axis:
        profile = epsf.data.sum(axis=a)
        itop = profile.argmax()
        poffset = profile[itop] / 2

        # Run a spline through the points, but subtract half of the peak value, and find the roots:
        # We have to use a cubic spline, since roots() is not supported for other splines
        # for some reason
        profile_intp = UnivariateSpline(np.arange(0, len(profile)), profile - poffset, k=3, s=0, ext=3)
        lr = profile_intp.roots()

        # Do some sanity checks on the ePSF:
        # It should pass 50% exactly twice and have the maximum inside that region.
        # I.e. it should be a single gaussian-like peak
        if len(lr) != 2 or itop < lr[0] or itop > lr[1]:
            logger.error(f"EPSF is not a single gaussian-like peak along axis {a}")
            epsf_ok = False
        else:
            axis_fwhm = lr[1] - lr[0]
            fwhms.append(axis_fwhm)

    return epsf_ok, fwhms
