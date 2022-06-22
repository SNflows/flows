import logging
from typing import Tuple, Any, Optional
from astropy.table import Table
import numpy as np
from bottleneck import nansum
from astropy.stats import sigma_clip
from astropy.modeling import models, fitting
import matplotlib.pyplot as plt

from .target import Target
from .zeropoint import sigma_from_Chauvenet, bootstrap_outlier
from .filters import get_reference_filter

logger = logging.getLogger(__name__)


def instrumental_mag(tab: Table, target: Target, make_fig: bool = False) -> Tuple[Table, Optional[plt.Figure],
                                                                                  Optional[plt.Axes]]:
    target_rows = tab['starid'] <= 0

    # Check that we got valid flux photometry:
    if np.any(~np.isfinite(tab[target_rows]['flux_psf'])) or np.any(~np.isfinite(tab[target_rows]['flux_psf_error'])):
        raise RuntimeError(f"Target:{target.name} flux is undefined.")

    # Convert PSF fluxes to magnitudes:
    mag_inst = -2.5 * np.log10(tab['flux_psf'])
    mag_inst_err = (2.5 / np.log(10)) * (tab['flux_psf_error'] / tab['flux_psf'])

    # Corresponding magnitudes in catalog:
    mag_catalog = tab[get_reference_filter(target.photfilter)]

    # Mask out things that should not be used in calibration:
    use_for_calibration = np.ones_like(mag_catalog, dtype='bool')
    use_for_calibration[target_rows] = False  # Do not use target for calibration
    use_for_calibration[~np.isfinite(mag_inst) | ~np.isfinite(mag_catalog)] = False


    # Just creating some short-hands:
    x = mag_catalog[use_for_calibration]
    y = mag_inst[use_for_calibration]
    yerr = mag_inst_err[use_for_calibration]
    weights = 1.0 / yerr ** 2

    if not any(use_for_calibration):
        raise RuntimeError("No calibration stars")

    # Fit linear function with fixed slope, using sigma-clipping:
    model = models.Linear1D(slope=1, fixed={'slope': True})
    fitter = fitting.FittingWithOutlierRemoval(fitting.LinearLSQFitter(), sigma_clip, sigma=3.0)
    best_fit, sigma_clipped = fitter(model, x, y, weights=weights)

    # Extract zero-point and estimate its error using a single weighted fit:
    # I don't know why there is not an error-estimate attached directly to the Parameter?
    zp = -1 * best_fit.intercept.value  # Negative, because that is the way zeropoints are usually defined

    weights[sigma_clipped] = 0  # Trick to make following expression simpler
    n_weights = len(weights.nonzero()[0])
    if n_weights > 1:
        zp_error = np.sqrt(n_weights * nansum(weights * (y - best_fit(x)) ** 2) / nansum(weights) / (n_weights - 1))
    else:
        zp_error = np.NaN
    logger.info('Leastsquare ZP = %.3f, ZP_error = %.3f', zp, zp_error)

    # Determine sigma clipping sigma according to Chauvenet method
    # But don't allow less than sigma = sigmamin, setting to 1.5 for now.
    # Should maybe be 2?
    sigmamin = 1.5
    sig_chauv = sigma_from_Chauvenet(len(x))
    sig_chauv = sig_chauv if sig_chauv >= sigmamin else sigmamin

    # Extract zero point and error using bootstrap method
    nboot = 1000
    logger.info('Running bootstrap with sigma = %.2f and n = %d', sig_chauv, nboot)
    pars = bootstrap_outlier(x, y, yerr, n=nboot, model=model, fitter=fitting.LinearLSQFitter, outlier=sigma_clip,
                             outlier_kwargs={'sigma': sig_chauv}, summary='median', error='bootstrap',
                             return_vals=False)

    zp_bs = pars['intercept'] * -1.0
    zp_error_bs = pars['intercept_error']

    logger.info('Bootstrapped ZP = %.3f, ZP_error = %.3f', zp_bs, zp_error_bs)

    # Check that difference is not large
    zp_diff = 0.4
    if np.abs(zp_bs - zp) >= zp_diff:
        logger.warning("Bootstrap and weighted LSQ ZPs differ by %.2f, "
                       "which is more than the allowed %.2f mag.", np.abs(zp_bs - zp), zp_diff)

    # Add calibrated magnitudes to the photometry table:
    tab['mag'] = mag_inst + zp_bs
    tab['mag_error'] = np.sqrt(mag_inst_err ** 2 + zp_error_bs ** 2)

    # Check that we got valid magnitude photometry:
    if not np.isfinite(tab[0]['mag']) or not np.isfinite(tab[0]['mag_error']):
        raise RuntimeError(f"Target:{target.name} magnitude is undefined.")

    # Update Meta-data:
    tab.meta['zp'] = zp_bs
    tab.meta['zp_error'] = zp_error_bs
    tab.meta['zp_diff'] = np.abs(zp_bs - zp)
    tab.meta['zp_error_weights'] = zp_error

    # Plot:
    if make_fig:
        mag_fig, mag_ax = plt.subplots(1, 1)
        mag_ax.errorbar(x, y, yerr=yerr, fmt='k.')
        mag_ax.scatter(x[sigma_clipped], y[sigma_clipped], marker='x', c='r')
        mag_ax.plot(x, best_fit(x), color='g', linewidth=3)
        mag_ax.set_xlabel('Catalog magnitude')
        mag_ax.set_ylabel('Instrumental magnitude')

        return tab, mag_fig, mag_ax
    return tab, None, None
