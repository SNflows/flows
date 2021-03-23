#import sys
import time

import numpy as np

#from scipy.ndimage import gaussian_filter

#from scipy.spatial import cKDTree
from scipy.interpolate import griddata

import photutils.psf

class EPSFBuilder(photutils.psf.EPSFBuilder):

    def _create_initial_epsf(self, stars):

        epsf = super()._create_initial_epsf(stars)
        epsf.origin = None

        X, Y = np.meshgrid(*map(np.arange, epsf.shape[::-1]))

        X = X / epsf.oversampling[0] - epsf.x_origin
        Y = Y / epsf.oversampling[1] - epsf.y_origin

        self._epsf_xy_grid = X, Y

        return epsf

    def _resample_residual(self, star, epsf):

        #max_dist = .5 / np.sqrt(np.sum(np.power(epsf.oversampling, 2)))

        #star_points = list(zip(star._xidx_centered, star._yidx_centered))
        #epsf_points = list(zip(*map(np.ravel, self._epsf_xy_grid)))

        #star_tree = cKDTree(star_points)
        #dd, ii = star_tree.query(epsf_points, distance_upper_bound=max_dist)
        #mask = np.isfinite(dd)

        #star_data = np.full_like(epsf.data, np.nan)
        #star_data.ravel()[mask] = star._data_values_normalized[ii[mask]]

        star_points = list(zip(star._xidx_centered, star._yidx_centered))
        star_data = griddata(star_points, star._data_values_normalized, self._epsf_xy_grid)

        return star_data - epsf._data

#    def _resample_residuals(self, stars, epsf):
#
#        residuals = super()._resample_residuals(stars, epsf)
#
#        import matplotlib.pyplot as plt
#
#        if epsf.data.any():
#
#            #i = np.isfinite(residuals).sum(axis=0, dtype=bool)
#            #data = epsf.data - gaussian_filter(epsf.data, sigma=1)
#            #factor = np.nanstd(data) / np.std(np.nanmedian(residuals[:,i], axis=0))
#            #residuals *= min(max(factor, 0.1), 1.0)
#
#            data = np.nanmedian(residuals, axis=0)
#            self._.set_data(data)
#            self._.set_clim([np.nanmin(data), np.nanmax(data)])
#            self._.cb.draw_all()
#
#        else:
#
#            data = np.nanmedian(residuals, axis=0)
#            self._ = plt.imshow(data)
#            self._.set_clim([np.nanmin(data), np.nanmax(data)])
#            self._.cb = plt.gcf().colorbar(self._)
#            plt.show(block=False)
#
#        plt.gcf().canvas.draw()
#        time.sleep(1)
#
#        return residuals

#    def _recenter_epsf(self, epsf, *args, **kwargs):
#
#        if not hasattr(self, 'dx_total'):
#
#            self.dx_total = []
#            self.dy_total = []
#
#        def profile(frame, event, arg):
#
#            global x, y, xcenter, ycenter
#
#            if event == "return" and frame.f_code.co_name == '_recenter_epsf':
#
#                x, xcenter = frame.f_locals['x'], frame.f_locals['xcenter']
#                y, ycenter = frame.f_locals['y'], frame.f_locals['ycenter']
#
#                self.dx_total.append(frame.f_locals['dx_total'])
#                self.dy_total.append(frame.f_locals['dy_total'])
#
#            return profile
#
#        sys.setprofile(profile)
#        super()._recenter_epsf(epsf, *args, **kwargs)
#        sys.setprofile(None)
#
#        dx_total = np.mean(self.dx_total[-2:], axis=0)
#        dy_total = np.mean(self.dy_total[-2:], axis=0)
#
#        epsf_data = epsf.evaluate(x=x, y=y, flux=1.0,
#                                  x_0=xcenter - .5*dx_total,
#                                  y_0=ycenter - .5*dy_total)
#
#        return epsf_data

    def __call__(self, *args, **kwargs):

        #import matplotlib
        #_backend = matplotlib.get_backend()
        #matplotlib.pyplot.switch_backend('TkAgg')

        t0 = time.time()

        epsf, stars = super().__call__(*args, **kwargs)

        epsf.fit_info = dict(
            n_iter = len(self._epsf),
            max_iters = self.maxiters,
            time = time.time() - t0,
        )

        #matplotlib.pyplot.switch_backend(_backend)

        return epsf, stars
