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
from .load_image import load_image, FlowsImage  # noqa: E402
from .coordinatematch import CoordinateMatch, WCS2  # noqa: E402
from .epsfbuilder import FlowsEPSFBuilder, verify_epsf  # noqa: E402
from .fileio import DirectoryProtocol, Directories  # noqa: E402
from .filters import get_reference_filter  # noqa: E402
from .target import Target  # noqa: E402

# @TODO: refactor load_image to separate modules.
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
    if diffimage_df:
        diffimage_path = diffimage_df.get('path', None)
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
        psfphot_tbl.insert_row(0, diff_psfphot_tbl[0])
        apphot_tbl.insert_row(0, diff_apphot_tbl[0])


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
        if diffimage_df is not None:
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


# def photometry(fileid, output_folder=None, attempt_imagematch=True, keep_diff_fixed=False, cm_timeout=None):
#     """
#     Run photometry.
#
#     Parameters:
#         fileid (int): File ID to process.
#         output_folder (str, optional): Path to directory where output should be placed.
#         attempt_imagematch (bool, optional): If no subtracted image is available, but a
#             template image is, should we attempt to run ImageMatch using standard settings.
#             Default=True.
#         keep_diff_fixed (bool, optional): Allow psf photometry to recenter when
#             calculating the flux for the difference image. Setting to True can help if diff
#             image has non-source flux in the region around the SN.
#         cm_timeout (float, optional): Timeout in seconds for the :class:`CoordinateMatch` algorithm.
#
#     .. codeauthor:: Rasmus Handberg <rasmush@phys.au.dk>
#     .. codeauthor:: Emir Karamehmetoglu <emir.k@phys.au.dk>
#     .. codeauthor:: Simon Holmbo <sholmbo@phys.au.dk>
#     """
#
#     # Settings:
#     ref_target_dist_limit = 10 * u.arcsec  # Reference star must be further than this away to be included
#
#     tic = default_timer()
#
#     # Use local copy of archive if configured to do so:
#     config = load_config()
#
#     # Get datafile dict from API:
#     datafile = api.get_datafile(fileid)
#     logger.debug("Datafile: %s", datafile)
#     targetid = datafile['targetid']
#     target_name = datafile['target_name']
#     photfilter = datafile['photfilter']
#
#     archive_local = config.get('photometry', 'archive_local', fallback=None)
#     if archive_local is not None:
#         datafile['archive_path'] = archive_local
#     if not os.path.isdir(datafile['archive_path']):
#         raise FileNotFoundError("ARCHIVE is not available: " + datafile['archive_path'])
#
#     # Get the catalog containing the target and reference stars:
#     # TODO: Include proper-motion to the time of observation
#     catalog = api.get_catalog(targetid, output='table')
#     target = catalog['target'][0]
#     target_coord = coords.SkyCoord(ra=target['ra'], dec=target['decl'], unit='deg', frame='icrs')
#
#     # Folder to save output:
#     if output_folder is None:
#         output_folder_root = config.get('photometry', 'output', fallback='.')
#         output_folder = os.path.join(output_folder_root, target_name, f'{fileid:05d}')
#     logger.info("Placing output in '%s'", output_folder)
#     os.makedirs(output_folder, exist_ok=True)
#
#     # The paths to the science image:
#     filepath = os.path.join(datafile['archive_path'], datafile['path'])
#
#     # TODO: Download datafile using API to local drive:
#     # TODO: Is this a security concern?
#     # if archive_local:
#     #	 api.download_datafile(datafile, archive_local)
#
#     # Translate photometric filter into table column:
#     ref_filter = get_reference_filter(photfilter)  # Fallback to 'g' if photfilter is not in FILTERS
#
#     # Load the image from the FITS file:
#     logger.info("Load image '%s'", filepath)
#     image = load_image(filepath, target_coord=target_coord)
#
#     references = catalog['references']
#     references.sort(ref_filter)
#
#     # Check that there actually are reference stars in that filter:
#     if allnan(references[ref_filter]):
#         raise ValueError("No reference stars found in current photfilter.")
#
#     # ==============================================================================================
#     # BARYCENTRIC CORRECTION OF TIME
#     # ==============================================================================================
#
#     ltt_bary = image.obstime.light_travel_time(target_coord, ephemeris='jpl')
#     image.obstime = image.obstime.tdb + ltt_bary
#
#     # ==============================================================================================
#     # BACKGROUND ESTIMATION
#     # ==============================================================================================
#
#     fig, ax = plt.subplots(1, 2, figsize=(20, 18))
#     plot_image(image.clean, ax=ax[0], scale='log', cbar='right', title='Image')
#     plot_image(image.mask, ax=ax[1], scale='linear', cbar='right', title='Mask')
#     fig.savefig(os.path.join(output_folder, 'original.png'), bbox_inches='tight')
#     plt.close(fig)
#
#     # Estimate image background:
#     # Not using image.clean here, since we are redefining the mask anyway
#     background = Background2D(image.clean, (128, 128), filter_size=(5, 5), sigma_clip=SigmaClip(sigma=3.0),
#                               bkg_estimator=SExtractorBackground(), exclude_percentile=50.0)
#
#     # Create background-subtracted image:
#     image.subclean = image.clean - background.background
#
#     # Plot background estimation:
#     fig, ax = plt.subplots(1, 3, figsize=(20, 6))
#     plot_image(image.clean, ax=ax[0], scale='log', title='Original')
#     plot_image(background.background, ax=ax[1], scale='log', title='Background')
#     plot_image(image.subclean, ax=ax[2], scale='log', title='Background subtracted')
#     fig.savefig(os.path.join(output_folder, 'background.png'), bbox_inches='tight')
#     plt.close(fig)
#
#     # TODO: Is this correct?!
#     image.error = calc_total_error(image.clean, background.background_rms, 1.0)
#
#     # Use sep to for soure extraction
#     sep_background = sep.Background(image.image, mask=image.mask)
#     objects = sep.extract(image.image - sep_background, thresh=5., err=sep_background.globalrms, mask=image.mask,
#                           deblend_cont=0.1, minarea=9, clean_param=2.0)
#
#     # Cleanup large arrays which are no longer needed:
#     del background, fig, ax, sep_background, ltt_bary
#     gc.collect()
#
#     # ==============================================================================================
#     # DETECTION OF STARS AND MATCHING WITH CATALOG
#     # ==============================================================================================
#
#     # Account for proper motion:
#     replace(references['pm_ra'], np.NaN, 0)
#     replace(references['pm_dec'], np.NaN, 0)
#     refs_coord = coords.SkyCoord(ra=references['ra'], dec=references['decl'], pm_ra_cosdec=references['pm_ra'],
#                                  pm_dec=references['pm_dec'], unit='deg', frame='icrs',
#                                  obstime=Time(2015.5, format='decimalyear'))
#
#     with warnings.catch_warnings():
#         warnings.simplefilter("ignore", ErfaWarning)
#         refs_coord = refs_coord.apply_space_motion(new_obstime=image.obstime)
#
#     # TODO: These need to be based on the instrument!
#     radius = 10
#     fwhm_guess = 6.0
#     fwhm_min = 3.5
#     fwhm_max = 18.0
#
#     # Clean extracted stars
#     masked_sep_xy, sep_mask, masked_sep_rsqs = refclean.force_reject_g2d(objects['x'], objects['y'], image,
#                                                                         radius=radius,
#                                                                          fwhm_guess=fwhm_guess, rsq_min=0.3,
#                                                                          fwhm_max=fwhm_max, fwhm_min=fwhm_min)
#
#     logger.info("Finding new WCS solution...")
#     head_wcs = str(WCS2.from_astropy_wcs(image.wcs))
#     logger.debug('Head WCS: %s', head_wcs)
#
#     # Solve for new WCS
#     cm = CoordinateMatch(xy=list(masked_sep_xy[sep_mask]), rd=list(zip(refs_coord.ra.deg, refs_coord.dec.deg)),
#                          xy_order=np.argsort(
#                              np.power(masked_sep_xy[sep_mask] - np.array(image.shape[::-1]) / 2, 2).sum(axis=1)),
#                          rd_order=np.argsort(target_coord.separation(refs_coord)), xy_nmax=100, rd_nmax=100,
#                          maximum_angle_distance=0.002)
#
#     # Set timeout par to infinity unless specified.
#     if cm_timeout is None:
#         cm_timeout = float('inf')
#     try:
#         i_xy, i_rd = map(np.array, zip(*cm(5, 1.5, timeout=cm_timeout)))
#     except TimeoutError:
#         logger.warning('TimeoutError: No new WCS solution found')
#     except StopIteration:
#         logger.warning('StopIterationError: No new WCS solution found')
#     else:
#         logger.info('Found new WCS')
#         image.wcs = fit_wcs_from_points(np.array(list(zip(*cm.xy[i_xy]))),
#                                         coords.SkyCoord(*map(list, zip(*cm.rd[i_rd])), unit='deg'))
#         del i_xy, i_rd
#
#     used_wcs = str(WCS2.from_astropy_wcs(image.wcs))
#     logger.debug('Used WCS: %s', used_wcs)
#
#     # Calculate pixel-coordinates of references:
#     xy = image.wcs.all_world2pix(list(zip(refs_coord.ra.deg, refs_coord.dec.deg)), 0)
#     references['pixel_column'], references['pixel_row'] = x, y = list(map(np.array, zip(*xy)))
#
#     # Clean out the references:
#     hsize = 10
#     ref_target_dist_limit = 10 * u.arcsec  # Reference star must be further than this away to be included
#     clean_references = references[(target_coord.separation(refs_coord) > ref_target_dist_limit) & (x > hsize) & (
#         x < (image.shape[1] - 1 - hsize)) & (y > hsize) & (y < (image.shape[0] - 1 - hsize))]
#
#     if not clean_references:
#         raise RuntimeError('No clean references in field')
#
#     # Calculate the targets position in the image:
#     target_pixel_pos = image.wcs.all_world2pix([(target['ra'], target['decl'])], 0)[0]
#
#     # Clean reference star locations
#     masked_fwhms, masked_ref_xys, rsq_mask, masked_rsqs = refclean.force_reject_g2d(clean_references['pixel_column'],
#                                                                                     clean_references['pixel_row'],
#                                                                                     image, get_fwhm=True, radius=radius,
#                                                                                     fwhm_guess=fwhm_guess,
#                                                                                     fwhm_max=fwhm_max,
#                                                                                     fwhm_min=fwhm_min, rsq_min=0.15)
#
#     # Use R^2 to more robustly determine initial FWHM guess.
#     # This cleaning is good when we have FEW references.
#     fwhm, fwhm_clean_references = refclean.clean_with_rsq_and_get_fwhm(masked_fwhms, masked_rsqs, clean_references,
#                                                                        min_fwhm_references=2, min_references=6,
#                                                                        rsq_min=0.15)
#     logger.info('Initial FWHM guess is %f pixels', fwhm)
#
#     # Create plot of target and reference star positions from 2D Gaussian fits.
#     fig, ax = plt.subplots(1, 1, figsize=(20, 18))
#     plot_image(image.subclean, ax=ax, scale='log', cbar='right', title=target_name)
#     ax.scatter(fwhm_clean_references['pixel_column'], fwhm_clean_references['pixel_row'], c='r', marker='o', alpha=0.3)
#     ax.scatter(masked_sep_xy[:, 0], masked_sep_xy[:, 1], marker='s', alpha=1.0, edgecolors='green', facecolors='none')
#     ax.scatter(target_pixel_pos[0], target_pixel_pos[1], marker='+', s=20, c='r')
#     fig.savefig(os.path.join(output_folder, 'positions_g2d.png'), bbox_inches='tight')
#     plt.close(fig)
#
#     # Final clean of wcs corrected references
#     logger.info("Number of references before final cleaning: %d", len(clean_references))
#     logger.debug('Masked R^2 values: %s', masked_rsqs[rsq_mask])
#     references = refclean.get_clean_references(clean_references, masked_rsqs, rsq_ideal=0.8)
#     logger.info("Number of references after final cleaning: %d", len(references))
#
#     # Create plot of target and reference star positions:
#     fig, ax = plt.subplots(1, 1, figsize=(20, 18))
#     plot_image(image.subclean, ax=ax, scale='log', cbar='right', title=target_name)
#     ax.scatter(references['pixel_column'], references['pixel_row'], c='r', marker='o', alpha=0.6)
#     ax.scatter(masked_sep_xy[:, 0], masked_sep_xy[:, 1], marker='s', alpha=0.6, edgecolors='green', facecolors='none')
#     ax.scatter(target_pixel_pos[0], target_pixel_pos[1], marker='+', s=20, c='r')
#     fig.savefig(os.path.join(output_folder, 'positions.png'), bbox_inches='tight')
#     plt.close(fig)
#
#     # Cleanup large arrays which are no longer needed:
#     del fig, ax, cm
#     gc.collect()
#
#     # ==============================================================================================
#     # CREATE EFFECTIVE PSF MODEL
#     # ==============================================================================================
#
#     # Make cutouts of stars using extract_stars:
#     # Scales with FWHM
#     size = int(np.round(29 * fwhm / 6))
#     size += 0 if size % 2 else 1  # Make sure it's a uneven number
#     size = max(size, 15)  # Never go below 15 pixels
#
#     # Extract stars sub-images:
#     xy = [tuple(masked_ref_xys[clean_references['starid'] == ref['starid']].data[0]) for ref in references]
#     with warnings.catch_warnings():
#         warnings.simplefilter('ignore', AstropyUserWarning)
#         stars = extract_stars(NDData(data=image.subclean.data, mask=image.mask), Table(np.array(xy), names=('x', 'y')),
#                               size=size + 6  # +6 for edge buffer
#                               )
#
#     logger.info("Number of stars input to ePSF builder: %d", len(stars))
#
#     # Plot the stars being used for ePSF:
#     imgnr = 0
#     nrows, ncols = 5, 5
#     for k in range(int(np.ceil(len(stars) / (nrows * ncols)))):
#         fig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=(20, 20), squeeze=True)
#         ax = ax.ravel()
#         for i in range(nrows * ncols):
#             if imgnr > len(stars) - 1:
#                 ax[i].axis('off')
#             else:
#                 plot_image(stars[imgnr], ax=ax[i], scale='log', cmap='viridis')  # FIXME (no x-ticks)
#             imgnr += 1
#
#         fig.savefig(os.path.join(output_folder, f'epsf_stars{k + 1:02d}.png'), bbox_inches='tight')
#         plt.close(fig)
#
#     # Build the ePSF:
#     epsf, stars = FlowsEPSFBuilder(oversampling=1, shape=1 * size,
#                                    fitter=EPSFFitter(fit_boxsize=max(int(np.round(1.5 * fwhm)), 5)),
#                                    recentering_boxsize=max(int(np.round(2 * fwhm)), 5), norm_radius=max(fwhm, 5),
#                                    maxiters=100, progress_bar=logger.isEnabledFor(logging.INFO))(stars)
#     logger.info('Built PSF model (%(n_iter)d/%(max_iters)d) in %(time).1f seconds', epsf.fit_info)
#
#     # Store which stars were used in ePSF in the table:
#     references['used_for_epsf'] = False
#     references['used_for_epsf'][[star.id_label - 1 for star in stars.all_good_stars]] = True
#     logger.info("Number of stars used for ePSF: %d", np.sum(references['used_for_epsf']))
#
#     fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 15))
#     plot_image(epsf.data, ax=ax1, cmap='viridis')
#
#     fwhms = []
#     bad_epsf_detected = False
#     for a, ax in ((0, ax3), (1, ax2)):
#         # Collapse the PDF along this axis:
#         profile = epsf.data.sum(axis=a)
#         itop = profile.argmax()
#         poffset = profile[itop] / 2
#
#         # Run a spline through the points, but subtract half of the peak value, and find the roots:
#         # We have to use a cubic spline, since roots() is not supported for other splines
#         # for some reason
#         profile_intp = UnivariateSpline(np.arange(0, len(profile)), profile - poffset, k=3, s=0, ext=3)
#         lr = profile_intp.roots()
#
#         # Plot the profile and spline:
#         x_fine = np.linspace(-0.5, len(profile) - 0.5, 500)
#         ax.plot(profile, 'k.-')
#         ax.plot(x_fine, profile_intp(x_fine) + poffset, 'g-')
#         ax.axvline(itop)
#         ax.set_xlim(-0.5, len(profile) - 0.5)
#
#         # Do some sanity checks on the ePSF:
#         # It should pass 50% exactly twice and have the maximum inside that region.
#         # I.e. it should be a single gaussian-like peak
#         if len(lr) != 2 or itop < lr[0] or itop > lr[1]:
#             logger.error("Bad PSF along axis %d", a)
#             bad_epsf_detected = True
#         else:
#             axis_fwhm = lr[1] - lr[0]
#             fwhms.append(axis_fwhm)
#             ax.axvspan(lr[0], lr[1], facecolor='g', alpha=0.2)
#
#     # Save the ePSF figure:
#     ax4.axis('off')
#     fig.savefig(os.path.join(output_folder, 'epsf.png'), bbox_inches='tight')
#     plt.close(fig)
#
#     # There was a problem with the ePSF:
#     if bad_epsf_detected:
#         raise RuntimeError("Bad ePSF detected.")
#
#     # Let's make the final FWHM the largest one we found:
#     fwhm = np.max(fwhms)
#     logger.info("Final FWHM based on ePSF: %f", fwhm)
#
#     # Cleanup large arrays which are no longer needed:
#     del fig, ax, stars, fwhms, profile_intp
#     gc.collect()
#
#     # ==============================================================================================
#     # COORDINATES TO DO PHOTOMETRY AT
#     # ==============================================================================================
#
#     coordinates = np.array([[ref['pixel_column'], ref['pixel_row']] for ref in references])
#
#     # Add the main target position as the first entry for doing photometry directly in the
#     # science image:
#     coordinates = np.concatenate(([target_pixel_pos], coordinates), axis=0)
#
#     # ==============================================================================================
#     # APERTURE PHOTOMETRY
#     # ==============================================================================================
#
#     # Define apertures for aperture photometry:
#     apertures = CircularAperture(coordinates, r=fwhm)
#     annuli = CircularAnnulus(coordinates, r_in=1.5 * fwhm, r_out=2.5 * fwhm)
#
#     apphot_tbl = aperture_photometry(image.subclean, [apertures, annuli], mask=image.mask, error=image.error)
#
#     logger.info('Aperture Photometry Success')
#     logger.debug("Aperture Photometry Table:\n%s", apphot_tbl)
#
#     # ==============================================================================================
#     # PSF PHOTOMETRY
#     # ==============================================================================================
#
#     # Create photometry object:
#     photometry_obj = BasicPSFPhotometry(group_maker=DAOGroup(fwhm), bkg_estimator=MedianBackground(), psf_model=epsf,
#                                         fitter=fitting.LevMarLSQFitter(), fitshape=size, aperture_radius=fwhm)
#
#     psfphot_tbl = photometry_obj(image=image.subclean, init_guesses=Table(coordinates, names=['x_0', 'y_0']))
#
#     logger.info('PSF Photometry Success')
#     logger.debug("PSF Photometry Table:\n%s", psfphot_tbl)
#
#     # ==============================================================================================
#     # TEMPLATE SUBTRACTION AND TARGET PHOTOMETRY
#     # ==============================================================================================
#
#     # Find the pixel-scale of the science image:
#     pixel_area = proj_plane_pixel_area(image.wcs.celestial)
#     pixel_scale = np.sqrt(pixel_area) * 3600  # arcsec/pixel
#     # print(image.wcs.celestial.cunit) % Doesn't work?
#     logger.info("Science image pixel scale: %f", pixel_scale)
#
#     diffimage = None
#     if datafile.get('diffimg') is not None:
#         diffimg_path = os.path.join(datafile['archive_path'], datafile['diffimg']['path'])
#         diffimg = load_image(diffimg_path)
#         diffimage = diffimg.image
#
#     elif attempt_imagematch and datafile.get('template') is not None:
#         # Run the template subtraction, and get back
#         # the science image where the template has been subtracted:
#         diffimage = run_imagematch(datafile, target, star_coord=coordinates, fwhm=fwhm, pixel_scale=pixel_scale)
#
#     # We have a diff image, so let's do photometry of the target using this:
#     if diffimage is not None:
#         # Include mask from original image:
#         diffimage = np.ma.masked_array(diffimage, image.mask)
#
#         # Create apertures around the target:
#         apertures = CircularAperture(target_pixel_pos, r=fwhm)
#         annuli = CircularAnnulus(target_pixel_pos, r_in=1.5 * fwhm, r_out=2.5 * fwhm)
#
#         # Create two plots of the difference image:
#         fig, ax = plt.subplots(1, 1, squeeze=True, figsize=(20, 20))
#         plot_image(diffimage, ax=ax, cbar='right', title=target_name)
#         ax.plot(target_pixel_pos[0], target_pixel_pos[1], marker='+', markersize=20, color='r')
#         fig.savefig(os.path.join(output_folder, 'diffimg.png'), bbox_inches='tight')
#         apertures.plot(axes=ax, color='r', lw=2)
#         annuli.plot(axes=ax, color='r', lw=2)
#         ax.set_xlim(target_pixel_pos[0] - 50, target_pixel_pos[0] + 50)
#         ax.set_ylim(target_pixel_pos[1] - 50, target_pixel_pos[1] + 50)
#         fig.savefig(os.path.join(output_folder, 'diffimg_zoom.png'), bbox_inches='tight')
#         plt.close(fig)
#
#         # Run aperture photometry on subtracted image:
#         target_apphot_tbl = aperture_photometry(diffimage, [apertures, annuli], mask=image.mask, error=image.error)
#
#         # Make target only photometry object if keep_diff_fixed = True
#         if keep_diff_fixed:
#             epsf.fixed.update({'x_0': True, 'y_0': True})
#
#             # TODO: Try iteraratively subtracted photometry
#             # Create photometry object:
#             photometry_obj = BasicPSFPhotometry(group_maker=DAOGroup(0.0001), bkg_estimator=MedianBackground(),
#                                                 psf_model=epsf, fitter=fitting.LevMarLSQFitter(), fitshape=size,
#                                                 aperture_radius=fwhm)
#
#         # Run PSF photometry on template subtracted image:
#         target_psfphot_tbl = photometry_obj(diffimage, init_guesses=Table(target_pixel_pos, names=['x_0', 'y_0']))
#
#         # Need to adjust table columns if x_0 and y_0 were fixed
#         if keep_diff_fixed:
#             target_psfphot_tbl['x_0_unc'] = 0.0
#             target_psfphot_tbl['y_0_unc'] = 0.0
#
#         # Combine the output tables from the target and the reference stars into one:
#         apphot_tbl = vstack([target_apphot_tbl, apphot_tbl], join_type='exact')
#         psfphot_tbl = vstack([target_psfphot_tbl, psfphot_tbl], join_type='exact')
#
#     # Build results table:
#     tab = references.copy()
#
#     row = {'starid': 0, 'ra': target['ra'], 'decl': target['decl'], 'pixel_column': target_pixel_pos[0],
#            'pixel_row': target_pixel_pos[1], 'used_for_epsf': False}
#     row.update([(k, np.NaN) for k in set(tab.keys()) - set(row) - {'gaia_variability'}])
#     tab.insert_row(0, row)
#
#     if diffimage is not None:
#         row['starid'] = -1
#         tab.insert_row(0, row)
#
#     indx_main_target = tab['starid'] <= 0
#
#     # Subtract background estimated from annuli:
#     flux_aperture = apphot_tbl['aperture_sum_0'] - (apphot_tbl['aperture_sum_1'] / annuli.area) * apertures.area
#     flux_aperture_error = np.sqrt(
#         apphot_tbl['aperture_sum_err_0'] ** 2 + (apphot_tbl['aperture_sum_err_1'] / annuli.area * apertures.area) ** 2)
#
#     # Add table columns with results:
#     tab['flux_aperture'] = flux_aperture / image.exptime
#     tab['flux_aperture_error'] = flux_aperture_error / image.exptime
#     tab['flux_psf'] = psfphot_tbl['flux_fit'] / image.exptime
#     tab['flux_psf_error'] = psfphot_tbl['flux_unc'] / image.exptime
#     tab['pixel_column_psf_fit'] = psfphot_tbl['x_fit']
#     tab['pixel_row_psf_fit'] = psfphot_tbl['y_fit']
#     tab['pixel_column_psf_fit_error'] = psfphot_tbl['x_0_unc']
#     tab['pixel_row_psf_fit_error'] = psfphot_tbl['y_0_unc']
#
#     # Check that we got valid photometry:
#     if np.any(~np.isfinite(tab[indx_main_target]['flux_psf'])) or np.any(
#         ~np.isfinite(tab[indx_main_target]['flux_psf_error'])):
#         raise RuntimeError("Target magnitude is undefined.")
#
#     # ==============================================================================================
#     # CALIBRATE
#     # ==============================================================================================
#
#     # Convert PSF fluxes to magnitudes:
#     mag_inst = -2.5 * np.log10(tab['flux_psf'])
#     mag_inst_err = (2.5 / np.log(10)) * (tab['flux_psf_error'] / tab['flux_psf'])
#
#     # Corresponding magnitudes in catalog:
#     mag_catalog = tab[ref_filter]
#
#     # Mask out things that should not be used in calibration:
#     use_for_calibration = np.ones_like(mag_catalog, dtype='bool')
#     use_for_calibration[indx_main_target] = False  # Do not use target for calibration
#     use_for_calibration[~np.isfinite(mag_inst) | ~np.isfinite(mag_catalog)] = False
#
#     # Just creating some short-hands:
#     x = mag_catalog[use_for_calibration]
#     y = mag_inst[use_for_calibration]
#     yerr = mag_inst_err[use_for_calibration]
#     weights = 1.0 / yerr ** 2
#
#     if not any(use_for_calibration):
#         raise RuntimeError("No calibration stars")
#
#     # Fit linear function with fixed slope, using sigma-clipping:
#     model = models.Linear1D(slope=1, fixed={'slope': True})
#     fitter = fitting.FittingWithOutlierRemoval(fitting.LinearLSQFitter(), sigma_clip, sigma=3.0)
#     best_fit, sigma_clipped = fitter(model, x, y, weights=weights)
#
#     # Extract zero-point and estimate its error using a single weighted fit:
#     # I don't know why there is not an error-estimate attached directly to the Parameter?
#     zp = -1 * best_fit.intercept.value  # Negative, because that is the way zeropoints are usually defined
#
#     weights[sigma_clipped] = 0  # Trick to make following expression simpler
#     n_weights = len(weights.nonzero()[0])
#     if n_weights > 1:
#         zp_error = np.sqrt(n_weights * nansum(weights * (y - best_fit(x)) ** 2) / nansum(weights) / (n_weights - 1))
#     else:
#         zp_error = np.NaN
#     logger.info('Leastsquare ZP = %.3f, ZP_error = %.3f', zp, zp_error)
#
#     # Determine sigma clipping sigma according to Chauvenet method
#     # But don't allow less than sigma = sigmamin, setting to 1.5 for now.
#     # Should maybe be 2?
#     sigmamin = 1.5
#     sig_chauv = sigma_from_Chauvenet(len(x))
#     sig_chauv = sig_chauv if sig_chauv >= sigmamin else sigmamin
#
#     # Extract zero point and error using bootstrap method
#     nboot = 1000
#     logger.info('Running bootstrap with sigma = %.2f and n = %d', sig_chauv, nboot)
#     pars = bootstrap_outlier(x, y, yerr, n=nboot, model=model, fitter=fitting.LinearLSQFitter, outlier=sigma_clip,
#                              outlier_kwargs={'sigma': sig_chauv}, summary='median', error='bootstrap',
#                              return_vals=False)
#
#     zp_bs = pars['intercept'] * -1.0
#     zp_error_bs = pars['intercept_error']
#
#     logger.info('Bootstrapped ZP = %.3f, ZP_error = %.3f', zp_bs, zp_error_bs)
#
#     # Check that difference is not large
#     zp_diff = 0.4
#     if np.abs(zp_bs - zp) >= zp_diff:
#         logger.warning("Bootstrap and weighted LSQ ZPs differ by %.2f, \
#             which is more than the allowed %.2f mag.", np.abs(zp_bs - zp), zp_diff)
#
#     # Add calibrated magnitudes to the photometry table:
#     tab['mag'] = mag_inst + zp_bs
#     tab['mag_error'] = np.sqrt(mag_inst_err ** 2 + zp_error_bs ** 2)
#
#     fig, ax = plt.subplots(1, 1)
#     ax.errorbar(x, y, yerr=yerr, fmt='k.')
#     ax.scatter(x[sigma_clipped], y[sigma_clipped], marker='x', c='r')
#     ax.plot(x, best_fit(x), color='g', linewidth=3)
#     ax.set_xlabel('Catalog magnitude')
#     ax.set_ylabel('Instrumental magnitude')
#     fig.savefig(os.path.join(output_folder, 'calibration.png'), bbox_inches='tight')
#     plt.close(fig)
#
#     # Check that we got valid photometry:
#     if not np.isfinite(tab[0]['mag']) or not np.isfinite(tab[0]['mag_error']):
#         raise RuntimeError("Target magnitude is undefined.")
#
#     # ==============================================================================================
#     # SAVE PHOTOMETRY
#     # ==============================================================================================
#
#     # Descriptions of columns:
#     tab['used_for_epsf'].description = 'Was object used for building ePSF?'
#     tab['mag'].description = 'Measured magnitude'
#     tab['mag'].unit = u.mag
#     tab['mag_error'].description = 'Error on measured magnitude'
#     tab['mag_error'].unit = u.mag
#     tab['flux_aperture'].description = 'Measured flux using aperture photometry'
#     tab['flux_aperture'].unit = u.count / u.second
#     tab['flux_aperture_error'].description = 'Error on measured flux using aperture photometry'
#     tab['flux_aperture_error'].unit = u.count / u.second
#     tab['flux_psf'].description = 'Measured flux using PSF photometry'
#     tab['flux_psf'].unit = u.count / u.second
#     tab['flux_psf_error'].description = 'Error on measured flux using PSF photometry'
#     tab['flux_psf_error'].unit = u.count / u.second
#     tab['pixel_column'].description = 'Location on image pixel columns'
#     tab['pixel_column'].unit = u.pixel
#     tab['pixel_row'].description = 'Location on image pixel rows'
#     tab['pixel_row'].unit = u.pixel
#     tab['pixel_column_psf_fit'].description = 'Measured location on image pixel columns from PSF photometry'
#     tab['pixel_column_psf_fit'].unit = u.pixel
#     tab[
#         'pixel_column_psf_fit_error'].description = 'Error on measured location on image pixel columns from PSF photometry'
#     tab['pixel_column_psf_fit_error'].unit = u.pixel
#     tab['pixel_row_psf_fit'].description = 'Measured location on image pixel rows from PSF photometry'
#     tab['pixel_row_psf_fit'].unit = u.pixel
#     tab['pixel_row_psf_fit_error'].description = 'Error on measured location on image pixel rows from PSF photometry'
#     tab['pixel_row_psf_fit_error'].unit = u.pixel
#
#     # Meta-data:
#     tab.meta['fileid'] = fileid
#     tab.meta['target_name'] = target_name
#     tab.meta['version'] = __version__
#     tab.meta['template'] = None if datafile.get('template') is None else datafile['template']['fileid']
#     tab.meta['diffimg'] = None if datafile.get('diffimg') is None else datafile['diffimg']['fileid']
#     tab.meta['photfilter'] = photfilter
#     tab.meta['fwhm'] = fwhm * u.pixel
#     tab.meta['pixel_scale'] = pixel_scale * u.arcsec / u.pixel
#     tab.meta['seeing'] = (fwhm * pixel_scale) * u.arcsec
#     tab.meta['obstime-bmjd'] = float(image.obstime.mjd)
#     tab.meta['zp'] = zp_bs
#     tab.meta['zp_error'] = zp_error_bs
#     tab.meta['zp_diff'] = np.abs(zp_bs - zp)
#     tab.meta['zp_error_weights'] = zp_error
#     tab.meta['head_wcs'] = head_wcs  # TODO: Are these really useful?
#     tab.meta['used_wcs'] = used_wcs  # TODO: Are these really useful?
#
#     # Filepath where to save photometry:
#     photometry_output = os.path.join(output_folder, 'photometry.ecsv')
#
#     # Write the final table to file:
#     tab.write(photometry_output, format='ascii.ecsv', delimiter=',', overwrite=True)
#
#     toc = default_timer()
#
#     logger.info("------------------------------------------------------")
#     logger.info("Success!")
#     logger.info("Main target: %f +/- %f", tab[0]['mag'], tab[0]['mag_error'])
#     logger.info("Photometry took: %.1f seconds", toc - tic)
#
#     return photometry_output
