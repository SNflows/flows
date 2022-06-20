#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Flows photometry code.

.. codeauthor:: Rasmus Handberg <rasmush@phys.au.dk>
.. codeauthor:: Emir Karamehmetoglu <emir.k@phys.au.dk>
.. codeauthor:: Simon Holmbo <sholmbo@phys.au.dk>
"""

import logging
import warnings
from timeit import default_timer
import numpy as np
from astropy.coordinates import SkyCoord
from bottleneck import allnan
from astropy.utils.exceptions import AstropyDeprecationWarning, AstropyUserWarning
import astropy.units as u
from astropy.stats import SigmaClip
from astropy.table import Table
from astropy.nddata import NDData
from astropy.modeling import fitting
from astropy.wcs.utils import proj_plane_pixel_area, fit_wcs_from_points

from tendrils import api
from tendrils.utils import load_config
from typing import Dict, List, Optional, Callable
from numpy.typing import ArrayLike

from .magnitudes import instrumental_mag
from .result_model import ResultsTable

warnings.simplefilter('ignore', category=AstropyDeprecationWarning)
from photutils import CircularAperture, CircularAnnulus, aperture_photometry  # noqa: E402
from photutils.psf import EPSFFitter, BasicPSFPhotometry, DAOGroup, extract_stars  # noqa: E402
from photutils.background import Background2D, SExtractorBackground, MedianBackground  # noqa: E402
from photutils.utils import calc_total_error  # noqa: E402

from . import reference_cleaning as refclean  # noqa: E402
from .plots import plt, plot_image  # noqa: E402
from .version import get_version  # noqa: E402
from .load_image import load_image  # noqa: E402
from .image import FlowsImage # noqa: E402
from .coordinatematch import CoordinateMatch, WCS2  # noqa: E402
from .epsfbuilder import FlowsEPSFBuilder, verify_epsf  # noqa: E402
from .fileio import DirectoryProtocol, Directories  # noqa: E402
from .filters import get_reference_filter  # noqa: E402
from .target import Target  # noqa: E402

__version__ = get_version(pep440=False)
logger = logging.getLogger(__name__)


def get_datafile(fileid: int) -> Dict:
    """
    Get datafile from API, log it, return.
    """
    # targetid = datafile['targetid']
    # target_name = datafile['target_name']
    # photfilter = datafile['photfilter']
    datafile = api.get_datafile(fileid)
    logger.debug("Datafile: %s", datafile)
    return datafile


def get_catalog(targetid: int) -> Dict:
    catalog = api.get_catalog(targetid, output='table')
    logger.debug(f"catalog obtained for target: {targetid}")
    return catalog


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


def correct_wcs(image: FlowsImage, references: refclean.References, target: Target,
                timeout: float = np.inf) -> FlowsImage:
    """
    Correct WCS of image to match the reference image.
    """
    # Start pre-cleaning
    sep_cleaner = refclean.ReferenceCleaner(image, references, rsq_min=0.3)

    # Use Source Extractor to make clean references
    sep_references_clean = sep_cleaner.make_sep_clean_references()

    # Find WCS
    logger.info("Finding new WCS solution...")
    head_wcs = str(WCS2.from_astropy_wcs(image.wcs))
    logger.debug('Head WCS: %s', head_wcs)
    # Solve for new WCS
    cm = CoordinateMatch(
        xy=list(sep_references_clean.masked),
        rd=list(zip(references.coords.ra.deg, references.coords.dec.deg)),
        xy_order=np.argsort(np.power(sep_references_clean.masked - np.array(image.shape[::-1]) / 2, 2).sum(axis=1)),
        rd_order=np.argsort(target.coords.separation(references.coords)),
        xy_nmax=100, rd_nmax=100, maximum_angle_distance=0.002)

    try:
        i_xy, i_rd = map(np.array, zip(*cm(5, 1.5, timeout=timeout)))
    except TimeoutError:
        logger.warning('TimeoutError: No new WCS solution found')
    except StopIteration:
        logger.warning('StopIterationError: No new WCS solution found')
    else:
        logger.info('Found new WCS')
        image.wcs = fit_wcs_from_points(np.array(list(zip(*cm.xy[i_xy]))),
                                        SkyCoord(*map(list, zip(*cm.rd[i_rd])), unit='deg'))
        del i_xy, i_rd

    logger.debug(f'Used WCS: {WCS2.from_astropy_wcs(image.wcs)}')
    return image

# @TODO: make photometry protocol
class Photometry:
    init_cutout_size: int = 29
    min_pixels: int = 15

    def __init__(self, image: FlowsImage, target: Target, fwhm_guess: float,
                 epsf_builder: FlowsEPSFBuilder = FlowsEPSFBuilder):
        self.image = image
        self.target = target
        self.fwhm = fwhm_guess
        self.epsf_builder = epsf_builder

    @property
    def star_size(self) -> int:
        # Make cutouts of stars using extract_stars:
        # Scales with FWHM
        size = int(np.round(self.init_cutout_size * self.fwhm / 6))
        size += 0 if size % 2 else 1  # Make sure it's odd
        size = max(size, self.min_pixels)  # Never go below 15 pixels
        return size

    def extract_star_cutouts(self, star_xys: np.ndarray) -> List[np.ma.MaskedArray]:
        """
        Extract star cutouts from the image.
        """
        # Extract stars from image
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', AstropyUserWarning)
            stars = extract_stars(NDData(data=self.image.subclean.data, mask=self.image.mask),
                                  Table(star_xys, names=('x', 'y')),
                                  size=self.star_size + 6  # +6 for edge buffer
                                  )
        logger.info("Number of stars input to ePSF builder: %d", len(stars))
        return stars

    def make_ePSF(self, stars):
        """
        Make an ePSF from the star cutouts.
        """
        # Build ePSF
        logger.info("Building ePSF...")
        builder = self.epsf_builder(
            oversampling=1, shape=1 * self.star_size,
            fitter=EPSFFitter(fit_boxsize=max(int(np.round(1.5 * self.fwhm)), 5)),
            recentering_boxsize=max(int(np.round(2 * self.fwhm)), 5),
            norm_radius=max(self.fwhm, 5), maxiters=100,
            progress_bar=logger.isEnabledFor(logging.INFO)
        )
        epsf, stars = builder(stars)

        logger.info(f"Built PSF model "
            f"{epsf.fit_info['n_iter']/epsf.fit_info['max_iters']} in {epsf.fit_info['time']} seconds")

        return epsf, stars


class PhotometryManager:
    """
    Implement a runner to shuffle data.
    """

    def __init__(self, target: Target,
                 directories: DirectoryProtocol):
        self.target = target
        self.directories = directories
        self.output_folder = directories.output_folder
        self.archive_local = directories.archive_local

    def load_science_image(self, image_path: str) -> FlowsImage:
        # The paths to the science image:
        science_image = self.directories.image_path(image_path)
        # Load the image from the FITS file:
        logger.info("Load image '%s'", science_image)
        return load_image(science_image, target_coord=self.target.coords)

    def get_filter(self):
        return get_reference_filter(self.target.photfilter)

    def load_references(self, catalog):
        use_filter = self.get_filter()
        references = catalog['references']
        references.sort(use_filter)
        # Check that there actually are reference stars in that filter:
        if allnan(references[use_filter]):
            raise ValueError("No reference stars found in current photfilter.")
        return refclean.References(table=references)


def calculate_appflux(apphot_tbl: Table, apertures: CircularAperture, annuli: CircularAnnulus) -> Table:
    """
    Calculate the aperture flux for the given apertures and annuli and append result to table.
    """
    # Subtract background estimated from annuli:
    bkg = (apphot_tbl['aperture_sum_1'] / annuli.area) * apertures.area
    apphot_tbl['flux_aperture'] = apphot_tbl['aperture_sum_0'] - bkg

    apphot_tbl['flux_aperture_error'] = np.sqrt(apphot_tbl['aperture_sum_err_0'] ** 2 +
                                                (apphot_tbl['aperture_sum_err_1'] / annuli.area * apertures.area) ** 2)
    return apphot_tbl


def apphot(coordinates: ArrayLike, image: FlowsImage, fwhm: float, use_raw: bool = False) -> Table:
    img = image.clean if use_raw else image.subclean
    apertures = CircularAperture(coordinates, r=fwhm)
    annuli = CircularAnnulus(coordinates, r_in=1.5 * fwhm, r_out=2.5 * fwhm)
    apphot_tbl = aperture_photometry(img, [apertures, annuli], mask=image.mask, error=image.error)
    return calculate_appflux(apphot_tbl, apertures, annuli)


def do_phot(fileid: int, cm_timeout: Optional[float] = None, make_plots: bool = True):
    # TODO: Timer should be moved out of this function.
    tic = default_timer()

    # Load config and set up directories
    config = load_config()
    directories = Directories(config)

    # Query API for datafile and catalogs
    datafile = get_datafile(fileid)
    catalog = get_catalog(datafile['targetid'])
    target = Target(catalog['target'][0]['ra'],
                    catalog['target'][0]['decl'],
                    name=datafile['target_name'],
                    id=datafile['targetid'],
                    photfilter=datafile['photfilter'])

    #  set output directories, creating if necessary
    directories.set_output_dirs(target.name, fileid)
    # science_image = directories.image_path(datafile['path'])

    # Set up photometry runner
    pr = PhotometryManager(target, directories)
    image = pr.load_science_image(datafile['path'])  # FlowsImage
    references = pr.load_references(catalog)  # Reference catalog
    references.make_sky_coords()  # Make sky coordinates
    references.propagate(image.obstime)  # get reference catalog at obstime

    # Estimate background and subtract from cleaned image
    bkg = FlowsBackground()
    image.subclean = bkg.background_subtract(image.clean)
    image.error = bkg.error(image.clean)

    # Correct WCS
    cm_timeout = cm_timeout if cm_timeout is not None else np.inf
    image = correct_wcs(image, references, target=target, timeout=cm_timeout)

    # Calculate pixel-coordinates of references:
    references.get_xy(image.wcs)
    references.make_pixel_columns()

    # Clean out the references:
    cleaner = refclean.ReferenceCleaner(image, references, rsq_min=0.15)
    # Reject references that are too close to target or edge of the image
    masked_references = cleaner.mask_edge_and_target(target.coords)
    if not masked_references.table:
        raise RuntimeError("No clean references in field")

    # Clean reference star locations
    clean_references, fwhm = cleaner.clean_references(masked_references) # Clean the masked references

    # EPSF creation
    phot = Photometry(image, target, fwhm)
    star_cutouts = phot.extract_star_cutouts(cleaner.gaussian_xys)
    epsf, stars = phot.make_ePSF(star_cutouts)
    epsf_ok, epsf_fwhms = verify_epsf(epsf)
    if not epsf_ok:
        raise RuntimeError("Bad ePSF detected.")

    # Store which stars were used in ePSF in the table:
    clean_references.table.add_column(col=[False], name='used_for_epsf')
    clean_references.table['used_for_epsf'][[star.id_label - 1 for star in stars.all_good_stars]] = True
    logger.info("Number of stars used for ePSF: %d", np.sum(clean_references.table['used_for_epsf']))

    fwhm = np.max(epsf_fwhms)  # Use the largest FWHM as new FWHM
    logger.info(f"Final FWHM based on ePSF: {fwhm}")

    # position in the image including target as row 0:
    target.calc_pixels(image.wcs)
    clean_references.add_target(target, starid=0)  # by default prepends target

    # Aperture photometry:
    apphot_tbl = apphot(clean_references.xy, image, fwhm)

    # PSF photometry:
    photometry_obj = BasicPSFPhotometry(group_maker=DAOGroup(fwhm), bkg_estimator=MedianBackground(), psf_model=epsf,
                                        fitter=fitting.LevMarLSQFitter(), fitshape=phot.star_size, aperture_radius=fwhm)
    psfphot_tbl = photometry_obj(image=image.subclean, init_guesses=Table(clean_references.xy, names=['x_0', 'y_0']))

    # Difference image photometry:
    diffimage_df = datafile.get('diffimg', None)
    diffimage = None
    if diffimage_df:
        diffimage_path = diffimage_df.get('path', None)
        logger.info("Found diffimg: %s, running difference photometry.", diffimage_path)
        if diffimage_path is None:
            logger.warning("No diffimage present but without path, skipping diffimage photometry")
        diffimage = load_image(directories.image_path(diffimage_path), target_coord=target.coords)
        diffimage.error = image.error
        diff_apphot_tbl = apphot(clean_references.xy[0], diffimage, fwhm, use_raw=True)
        diff_psfphot_tbl = photometry_obj(image=diffimage.clean, init_guesses=Table(clean_references.xy[0],
                                                                                    names=['x_0', 'y_0']))

        # Store the difference image photometry on row 0 of the table.
        # This pushes the un-subtracted target photometry to row 1.
        clean_references.add_target(target, starid=-1)
        psfphot_tbl.insert_row(0, dict(diff_psfphot_tbl[0]))
        apphot_tbl.insert_row(0, dict(diff_apphot_tbl[0]))


    # TODO: This should be moved to the photometry manager.
    # Build results table:
    results_table = ResultsTable.make_results_table(clean_references.table, apphot_tbl, psfphot_tbl, image)

    # Todo: refactor.
    # Get instrumental magnitude (currently we do too much here).
    results_table, (mag_fig, mag_ax) = instrumental_mag(results_table, target)

    # Add metadata to the results table:
    results_table.meta['fileid'] = fileid
    results_table.meta['target_name'] = target.name
    results_table.meta['version'] = __version__
    results_table.meta['template'] = None if datafile.get('template') is None else datafile['template']['fileid']
    results_table.meta['diffimg'] = None if datafile.get('diffimg') is None else datafile['diffimg']['fileid']
    results_table.meta['photfilter'] = target.photfilter
    results_table.meta['fwhm'] = fwhm * u.pixel
    # Find the pixel-scale of the science image
    pixel_area = proj_plane_pixel_area(image.wcs.celestial)
    pixel_scale = np.sqrt(pixel_area) * 3600  # arcsec/pixel
    logger.info("Science image pixel scale: %f", pixel_scale)
    # Keep adding metadata now that we have the pixel scale.
    results_table.meta['pixel_scale'] = pixel_scale * u.arcsec / u.pixel
    results_table.meta['seeing'] = (fwhm * pixel_scale) * u.arcsec
    results_table.meta['obstime-bmjd'] = float(image.obstime.mjd)
    results_table.meta['used_wcs'] = str(image.wcs)

    # Save the results table:
    # TODO: Photometry should RETURN a table, not save it.
    # TODO: Saving should be offloaded to the caller or to a parameter at least.
    results_table.write(directories.photometry_path, format='ascii.ecsv', delimiter=',', overwrite=True)

    # Log result and time taken:
    # TODO: This should be logged by the calling function.
    logger.info("------------------------------------------------------")
    logger.info("Success!")
    logger.info("Main target: %f +/- %f", results_table[0]['mag'], results_table[0]['mag_error'])
    logger.info("Photometry took: %.1f seconds", default_timer() - tic)

    # Plotting. TODO: refactor.
    # These are here for now for backwards compatibility.
    if make_plots:
        # Plot the image:
        fig, ax = plt.subplots(1, 2, figsize=(20, 18))
        plot_image(image.clean, ax=ax[0], scale='log', cbar='right', title='Image')
        plot_image(image.mask, ax=ax[1], scale='linear', cbar='right', title='Mask')
        fig.savefig(directories.save_as('original.png'), bbox_inches='tight')
        plt.close(fig)

        # Plot background estimation:
        fig, ax = plt.subplots(1, 3, figsize=(20, 6))
        plot_image(image.clean, ax=ax[0], scale='log', title='Original')
        plot_image(bkg.background, ax=ax[1], scale='log', title='Background')
        plot_image(image.subclean, ax=ax[2], scale='log', title='Background subtracted')
        fig.savefig(directories.save_as('background.png'), bbox_inches='tight')
        plt.close(fig)

        # Create plot of target and reference star positions:
        fig, ax = plt.subplots(1, 1, figsize=(20, 18))
        plot_image(image.subclean, ax=ax, scale='log', cbar='right', title=target.name)
        ax.scatter(references.table['pixel_column'], references.table['pixel_row'], c='r', marker='o', alpha=0.6)
        ax.scatter(cleaner.gaussian_xys[:, 0], cleaner.gaussian_xys[:, 1], marker='s', alpha=0.6, edgecolors='green',
                   facecolors='none')
        ax.scatter(target.pixel_column, target.pixel_row, marker='+', s=20, c='r')
        fig.savefig(directories.save_as('positions.png'), bbox_inches='tight')
        plt.close(fig)

        # Plot EPSF
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 15))
        plot_image(epsf.data, ax=ax1, cmap='viridis')
        for a, ax in ((0, ax3), (1, ax2)):
            profile = epsf.data.sum(axis=a)
            ax.plot(profile, 'k.-')
            ax.axvline(profile.argmax())
            ax.set_xlim(-0.5, len(profile) - 0.5)
        ax4.axis('off')
        fig.savefig(directories.save_as('epsf.png'), bbox_inches='tight')
        plt.close(fig)
        del ax, ax1, ax2, ax3, ax4

        # Create two plots of the difference image:
        if diffimage_df is not None and diffimage is not None:
            fig, ax = plt.subplots(1, 1, squeeze=True, figsize=(20, 20))
            plot_image(diffimage.clean, ax=ax, cbar='right', title=target.name)
            ax.plot(target.pixel_column, target.pixel_row, marker='+', markersize=20, color='r')
            fig.savefig(directories.save_as('diffimg.png'), bbox_inches='tight')
            #apertures.plot(axes=ax, color='r', lw=2)
            #annuli.plot(axes=ax, color='r', lw=2)
            ax.set_xlim(target.pixel_column - 50, target.pixel_column + 50)
            ax.set_ylim(target.pixel_row - 50, target.pixel_row + 50)
            fig.savefig(directories.save_as('diffimg_zoom.png'), bbox_inches='tight')
            plt.close(fig)

        # Calibration (handled in magnitudes.py).
        mag_fig.savefig(directories.save_as('calibration.png'), bbox_inches='tight')
        plt.close(mag_fig)

    # TODO: Return results table or photometry object or an event signature, not the path to the file..
    return directories.photometry_path
