#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Flows photometry code.

.. codeauthor:: Rasmus Handberg <rasmush@phys.au.dk>
.. codeauthor:: Emir Karamehmetoglu <emir.k@phys.au.dk>
.. codeauthor:: Simon Holmbo <sholmbo@phys.au.dk>
"""
import warnings
from copy import copy
from timeit import default_timer
from typing import Dict, List, Optional, Callable, Tuple, Union, TypeVar, Any

import numpy as np
from numpy.typing import ArrayLike, NDArray
from astropy.utils.exceptions import AstropyDeprecationWarning, AstropyUserWarning
import astropy.units as u
from astropy.table import Table
from astropy.nddata import NDData
from astropy.modeling import fitting
from astropy.wcs.utils import proj_plane_pixel_area
import multiprocessing
#from tendrils import api

from .magnitudes import instrumental_mag
from .result_model import ResultsTable

warnings.simplefilter('ignore', category=AstropyDeprecationWarning)
from photutils import CircularAperture, CircularAnnulus, aperture_photometry  # noqa: E402
from photutils.psf import EPSFFitter, BasicPSFPhotometry, DAOGroup, extract_stars  # noqa: E402
from photutils.background import MedianBackground  # noqa: E402
import photutils  # noqa: E402

from .reference_cleaning import References, ReferenceCleaner, InitGuess  # noqa: E402
from .plots import plt, plot_image  # noqa: E402
from .version import get_version  # noqa: E402
from .image import FlowsImage  # noqa: E402
from .coordinatematch import correct_wcs  # noqa: E402
from .epsfbuilder import FlowsEPSFBuilder, verify_epsf  # noqa: E402
from .fileio import DirectoryProtocol, IOManager  # noqa: E402
from .target import Target  # noqa: E402
from .background import FlowsBackground  # noqa: E402
from .utilities import create_logger  # noqa: E402

__version__ = get_version(pep440=False)

PhotutilsBackground = TypeVar('PhotutilsBackground', bound=photutils.background.core.BackgroundBase)
logger = create_logger()


# def get_datafile(fileid: int) -> Dict:
#     """
#     Get datafile from API, log it, return.
#     """
#     datafile = api.get_datafile(fileid)
#     logger.debug("Datafile: %s", datafile)
#     return datafile


# def get_catalog(targetid: int) -> Dict:
#     catalog = api.get_catalog(targetid, output='table')
#     logger.debug(f"catalog obtained for target: {targetid}")
#     return catalog


class PSFBuilder:
    init_cutout_size: int = 29
    min_pixels: int = 15

    def __init__(self, image: FlowsImage, target: Target, fwhm_guess: float,
                 epsf_builder: FlowsEPSFBuilder = FlowsEPSFBuilder):
        self.image = image
        self.target = target
        self.fwhm = fwhm_guess
        self.epsf_builder = epsf_builder

        # To be updated later.
        self.epsf = None

    # @TODO: Move to PhotometryMediator
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

    def make_epsf(self, stars):
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
            progress_bar=multiprocessing.parent_process() is None and logger.getEffectiveLevel() <= 20
        )
        epsf, stars = builder(stars)

        logger.info(f"Built PSF model "
                    f"{epsf.fit_info['n_iter'] / epsf.fit_info['max_iters']} in {epsf.fit_info['time']} seconds")

        return epsf, stars


class Photometry:

    def __init__(self, photometry_obj: Optional[BasicPSFPhotometry] = None):
        self.photometry_obj = photometry_obj

    @staticmethod
    def calculate_appflux(apphot_tbl: Table, apertures: CircularAperture, annuli: CircularAnnulus) -> Table:
        """
        Calculate the aperture flux for the given apertures and annuli and append result to table.
        """
        # Subtract background estimated from annuli:
        bkg = (apphot_tbl['aperture_sum_1'] / annuli.area) * apertures.area
        apphot_tbl['flux_aperture'] = apphot_tbl['aperture_sum_0'] - bkg

        apphot_tbl['flux_aperture_error'] = np.sqrt(apphot_tbl['aperture_sum_err_0'] ** 2 +
                                                    (apphot_tbl[
                                                         'aperture_sum_err_1'] / annuli.area * apertures.area) ** 2)
        return apphot_tbl

    def apphot(self, coordinates: ArrayLike, image: FlowsImage, fwhm: float, use_raw: bool = False) -> Table:
        img = image.clean if use_raw else image.subclean
        apertures = CircularAperture(coordinates, r=fwhm)
        annuli = CircularAnnulus(coordinates, r_in=1.5 * fwhm, r_out=2.5 * fwhm)
        apphot_tbl = aperture_photometry(img, [apertures, annuli], mask=image.mask, error=image.error)
        return self.calculate_appflux(apphot_tbl, apertures, annuli)

    def create_photometry_object(self, fwhm: Union[float, u.Quantity], psf_model: photutils.psf.EPSFModel,
                                 fitsize: Union[int, Tuple[int]], fitter: Callable = fitting.LevMarLSQFitter(),
                                 bkg: PhotutilsBackground = MedianBackground()):
        self.photometry_obj = BasicPSFPhotometry(group_maker=DAOGroup(fwhm), bkg_estimator=bkg, psf_model=psf_model,
                                                 fitter=fitter, fitshape=fitsize, aperture_radius=fwhm)

    def psfphot(self, image: ArrayLike, init_table: Table) -> Tuple[BasicPSFPhotometry, Table]:
        """PSF photometry on init guesses table/row.
        """
        if self.photometry_obj is None:
            raise ValueError('Photometry object not initialized.')
        # logger.info(f"{init_table}")
        output: Table = self.photometry_obj(image=image, init_guesses=init_table)
        return self.photometry_obj, output

    @staticmethod
    def rescale_flux_error(phot_tables: Dict[int, Table],
                           flux: float, flux_err: float, exptime: float) -> tuple[float, float]:
        """Rescale the error using input phot_tables dict with keys as input fit shapes and values as tables/rows
        returning the new error estimate and its fit shape."""
        select_first_row = isinstance(list(phot_tables.values())[0], Table)

        for fit_shape, row in phot_tables.items():
            row = row[0] if select_first_row else row
            new_err = row['flux_unc'] / exptime
            new_flux = row['flux_fit'] / exptime
            if new_flux <= flux + flux_err:
                return fit_shape, new_err
        logger.warning("Rescaled psf flux errors do not overlap input flux + error to 1 sigma, using original error.")
        return 0, flux_err

    @staticmethod
    def get_fit_shapes(fwhm: Union[float, int], star_size: int, fwhm_min: int = 2,
                       fwhm_max: int = 4) -> NDArray[np.int_]:
        if star_size / fwhm < fwhm_min:
            return np.atleast_1d(np.array(star_size))
        fit_shapes = np.arange(int(fwhm_min * fwhm) - 1, min(int(fwhm_max * fwhm), star_size), 1)
        return fit_shapes[fit_shapes % 2 == 1]  # odd


class PhotometryManager:

    def __init__(self, target: Target, image: FlowsImage, bkg: FlowsBackground,
                 references: References, directories: DirectoryProtocol,
                 fwhm: Optional[float] = None,
                 psf_builder: Optional[PSFBuilder] = None,
                 cleaner: Optional[ReferenceCleaner] = None,
                 diffimage: Optional[FlowsImage] = None):
        self.target = target
        self.image = image
        self.bkg = bkg
        self.references = references
        self.directories = directories
        # Initially possibly None:
        self.fwhm = fwhm
        self.psf_builder = psf_builder
        self.cleaner = cleaner
        self.diffimage = diffimage
        # To be updated later.
        self.clean_references = references
        self.clean_references_with_diff = references
        self.init_guesses = None
        self.init_guesses_diff = None
        self.photometry = Photometry()
        self.diff_im_exists = diffimage is not None
        self.results_table = None

    def propogate_references(self):
        self.references.make_sky_coords()  # Make sky coordinates
        self.references.propagate(self.image.obstime)  # get reference catalog at obstime

    def background_subtract(self):
        self.image.subclean = self.bkg.background_subtract(self.image.clean)
        self.image.error = self.bkg.error(self.image.clean)
        if self.diff_im_exists:
            self.diffimage.error = self.image.error

    def recalculate_image_wcs(self, cm_timeout: float):
        self.image = correct_wcs(self.image, self.references, target=self.target, timeout=cm_timeout)

    def calculate_pixel_coordinates(self):
        # Calculate pixel-coordinates of references:
        self.references.get_xy(self.image.wcs)
        self.references.make_pixel_columns()

    def clean_reference_stars(self, rsq_min: float = 0.15):
        # Clean out the references:
        self.cleaner = ReferenceCleaner(self.image, self.references, rsq_min=rsq_min)
        # Reject references that are too close to target or edge of the image
        masked_references = self.cleaner.mask_edge_and_target(self.target.coords)
        if not masked_references.table:
            raise RuntimeError("No clean references in field")

        # Clean masked reference star locations
        self.clean_references, self.fwhm = self.cleaner.clean_references(masked_references)

    def create_psf_builder(self) -> PSFBuilder:
        return PSFBuilder(self.image, self.target, self.fwhm)

    def update_reference_with_epsf(self, stars):
        # Store which stars were used in ePSF in the table:
        self.clean_references.table.add_column(col=[False], name='used_for_epsf')
        self.clean_references.table['used_for_epsf'][[star.id_label - 1 for star in stars.all_good_stars]] = True
        logger.info("Number of stars used for ePSF: %d", np.sum(self.clean_references.table['used_for_epsf']))

    def create_epsf(self, psf_builder: PSFBuilder = None):
        if psf_builder is None:
            psf_builder = self.create_psf_builder()

        # EPSF creation
        star_cutouts = psf_builder.extract_star_cutouts(self.cleaner.gaussian_xys)
        epsf, stars = psf_builder.make_epsf(star_cutouts)
        epsf_ok, epsf_fwhms = verify_epsf(epsf)
        if not epsf_ok:
            raise RuntimeError("Bad ePSF detected.")
        psf_builder.epsf = epsf

        # Use the largest FWHM as new FWHM
        fwhm = np.max(epsf_fwhms)
        logger.info(f"Final FWHM based on ePSF: {fwhm}")
        psf_builder.fwhm = fwhm

        # Update state
        self.fwhm = fwhm
        self.psf_builder = psf_builder
        self.update_reference_with_epsf(stars)

    def prepare_target_and_references_for_photometry(self):
        # position in the image including target as row 0:
        self.target.calc_pixels(self.image.wcs)
        self.clean_references.add_target(self.target, starid=0)  # by default prepends target
        self.init_guesses = InitGuess(self.clean_references, target_row=0)
        if self.diff_im_exists:
            self.clean_references_with_diff = copy(self.clean_references)
            self.clean_references_with_diff.add_target(self.target, starid=-1)
        self.init_guesses_diff = InitGuess(self.clean_references_with_diff, target_row=1, diff_row=0)

    def apphot(self) -> Table:
        # apphot_tbl.insert_row(0, dict(diff_apphot_tbl[0]))
        apphot_tbl = self.photometry.apphot(self.clean_references.xy, self.image, self.fwhm)
        if self.diff_im_exists:
            # Add diff image photometry:
            diff_tbl = self.photometry.apphot(self.clean_references.xy[0], self.diffimage, self.fwhm, use_raw=True)
            apphot_tbl.insert_row(0, dict(diff_tbl[0]))
        return apphot_tbl

    def psfphot(self, fit_shape: Optional[Union[int, Tuple[int, int]]] = None) -> Table:
        fit_shape = self.psf_builder.star_size if fit_shape is None else fit_shape
        # PSF photometry:
        self.photometry.create_photometry_object(fwhm=self.fwhm, psf_model=self.psf_builder.epsf, fitsize=fit_shape)
        psfphot_tbl = self.raw_psf_phot()
        if self.diff_im_exists:
            # Add diff image photometry:
            diff_tbl = self.diff_psf_phot()
            psfphot_tbl.insert_row(0, dict(diff_tbl[0]))
        return psfphot_tbl

    def raw_psf_phot(self, init_guess: Optional[Table] = None) -> Table:
        init_guess = self.init_guesses.init_guess_full if init_guess is None else init_guess
        _, psf_tbl = self.photometry.psfphot(image=self.image.subclean, init_table=init_guess)
        return psf_tbl

    def diff_psf_phot(self) -> Table:
        _, psf_tbl = self.photometry.psfphot(image=self.diffimage.clean,
                                             init_table=self.init_guesses_diff.init_guess_diff)
        return psf_tbl

    def rescale_uncertainty(self, psfphot_tbl: Table, dynamic: bool = True, 
                            static_fwhm: float = 2.5, epsilon_mag: float = 0.004,
                            ensure_greater: bool = True):
        """
        Rescale the uncertainty of the PSF photometry using a variable fitsize.

        Parameters
        ----------
        psfphot_tbl : Table
            Table of PSF photometry.
        dynamic : bool
            Dynamically decide FWHM multiple for rescaling.
        static_fwhm : float
            FWHM multiple to use incase dynamic fails or don't want to use it. Default 2.5 determined empirically.
        epsilon_mag : float
            Small magnitude change within which new and old uncertainties are considered the same. 
            Should be smaller than ~1/2 the expected uncertainty.
        """
        # Rescale psf errors from fit iteratively
        fit_shapes = self.photometry.get_fit_shapes(self.fwhm, self.psf_builder.star_size)
        fit_shape = int(static_fwhm * self.fwhm)
        fit_shape = fit_shape if fit_shape % 2 == 1 else fit_shape + 1
        if dynamic and len(fit_shapes) > 1:
            _phot_tables_dict = {}
            for fitshape in fit_shapes:
                self.photometry.create_photometry_object(
                    fwhm=self.fwhm, psf_model=self.psf_builder.epsf, fitsize=fitshape)
                if self.diff_im_exists:
                    _table = self.diff_psf_phot()
                _table = self.raw_psf_phot(self.init_guesses.init_guess_target)
                if "flux_unc" in _table.colnames:
                    _phot_tables_dict[fitshape] = _table

            if len(_phot_tables_dict) == 0:
                logger.warning("No PSF errors found for dynamic rescaling, trying static.")
                dynamic = False
            else:
                # Find the fit shape elbow:
                flux = psfphot_tbl[0]['flux_fit']
                flux_err = psfphot_tbl[0]['flux_unc']
                exptime = self.image.exptime
                dynamic_fit_shape, new_err = self.photometry.rescale_flux_error(_phot_tables_dict, flux, flux_err,
                                                                                exptime)
                fit_shape = dynamic_fit_shape if dynamic_fit_shape != 0 else fit_shape

        # Recalculate all reference uncertainties using new fitsize:
        logger.info(f"Recalculating all reference uncertainties using new fitsize:"
                    f" {fit_shape} pixels, ({fit_shape/self.fwhm if dynamic else static_fwhm :.2} * FWHM).")
        psfphot_tbl_rescaled = self.psfphot(fit_shape)
        if psfphot_tbl['flux_unc'][0] > psfphot_tbl_rescaled['flux_unc'][0] + epsilon_mag and ensure_greater:
            logger.info("Recalculated uncertainties were smaller than original and ``ensure_greater`` was True:"
                        "Not using rescaled uncertainties for the SN.")
            psfphot_tbl['flux_unc'][1:] = psfphot_tbl_rescaled['flux_unc'][1:]
            return psfphot_tbl

        psfphot_tbl['flux_unc'] = psfphot_tbl_rescaled['flux_unc']
        return psfphot_tbl

    def make_result_table(self, psfphot_tbl: Table, apphot_tbl: Table):
        # Build results table:
        clean_references = self.clean_references_with_diff if self.diff_im_exists else self.clean_references
        self.results_table = ResultsTable.make_results_table(clean_references.table, apphot_tbl, psfphot_tbl,
                                                             self.image)

    def calculate_mag(self, make_plot: bool = False) -> Tuple[Optional[plt.Figure], Optional[plt.Axes]]:
        if self.results_table is None:
            raise ValueError("Results table is not initialized. Run photometry first.")
        # Todo: refactor.
        # Get instrumental magnitude (currently we do too much here).
        results_table, mag_fig, mag_ax = instrumental_mag(self.results_table, self.target, make_plot)
        self.results_table = results_table
        return mag_fig, mag_ax

    def calculate_pixel_scale(self):
        # Find the pixel-scale of the science image
        pixel_area = proj_plane_pixel_area(self.image.wcs.celestial)
        pixel_scale = np.sqrt(pixel_area) * 3600  # arcsec/pixel
        logger.info("Science image pixel scale: %f", pixel_scale)
        return pixel_scale

    def add_metadata(self):
        # Add metadata to the results table:
        self.results_table.meta['fileid'] = self.image.fid
        self.results_table.meta['target_name'] = self.target.name
        self.results_table.meta['version'] = __version__
        self.results_table.meta['template'] = self.image.template_fid
        self.results_table.meta['diffimg'] = self.diffimage.fid if self.diff_im_exists else None
        self.results_table.meta['photfilter'] = self.target.photfilter
        self.results_table.meta['fwhm'] = self.fwhm * u.pixel
        pixel_scale = self.calculate_pixel_scale()
        self.results_table.meta['pixel_scale'] = pixel_scale * u.arcsec / u.pixel
        self.results_table.meta['seeing'] = (self.fwhm * pixel_scale) * u.arcsec
        self.results_table.meta['obstime-bmjd'] = float(self.image.obstime.mjd)
        self.results_table.meta['used_wcs'] = str(self.image.wcs)

    @classmethod
    def create_from_fid(cls, fid: int, directories: Optional[DirectoryProtocol] = None,
                        create_directories: bool = True, datafile: Optional[Dict] = None) -> 'PhotometryManager':
        """
        Create a Photometry object from a fileid.
        """
        io = IOManager.from_fid(fid, directories=directories, create_directories=create_directories, datafile=datafile)
        return PhotometryManager(target=io.target, image=io.load_science_image(), bkg=FlowsBackground(),
                                 references=io.load_references(), directories=io.directories,
                                 diffimage=io.load_diff_image())


def do_phot(fileid: int, cm_timeout: Optional[float] = None, make_plots: bool = True,
            directories: Optional[DirectoryProtocol] = None, datafile: Optional[Dict[str, Any]] = None,
            rescale: bool = True, rescale_dynamic: bool = True) -> ResultsTable:
    # Set up photometry runner
    pm = PhotometryManager.create_from_fid(fileid, directories=directories, datafile=datafile, create_directories=True)

    # Set up photometry
    pm.propogate_references()
    pm.background_subtract()
    pm.recalculate_image_wcs(cm_timeout if cm_timeout is not None else np.inf)
    pm.calculate_pixel_coordinates()
    pm.clean_reference_stars()
    pm.create_epsf()  # Using default EPSF constructor, could pass in custom.
    pm.prepare_target_and_references_for_photometry()  # Add target to reference table

    # Do photometry
    apphot_tbl = pm.apphot()
    # Verify uncertainty exists after PSF phot:
    psfphot_tbl = ResultsTable.verify_uncertainty_column(pm.psfphot())  
    if rescale:  # Rescale uncertainties
        psfphot_tbl = pm.rescale_uncertainty(psfphot_tbl, dynamic=rescale_dynamic) 

    # Build results table and calculate magnitudes
    pm.make_result_table(psfphot_tbl, apphot_tbl)
    mag_fig, mag_ax = pm.calculate_mag(make_plot=make_plots)
    pm.add_metadata()

    if make_plots:
        do_plots(pm, mag_fig, mag_ax)
    return pm.results_table


def timed_photometry(fileid: int, cm_timeout: Optional[float] = None, make_plots: bool = True,
                     directories: Optional[DirectoryProtocol] = None, save: bool = True,
                     datafile: Optional[Dict[str, Any]] = None, rescale: bool = True,
                     rescale_dynamic: bool = True) -> ResultsTable:
    # TODO: Timer should be moved out of this function.
    tic = default_timer()
    results_table = do_phot(fileid, cm_timeout, make_plots, directories, datafile, rescale, rescale_dynamic)

    # Save the results table:
    if save:
        results_table.write(directories.photometry_path, format='ascii.ecsv', delimiter=',', overwrite=True)

    # Log result and time taken:
    logger.info("------------------------------------------------------")
    logger.info("Success!")
    logger.info("Main target: %f +/- %f", results_table[0]['mag'], results_table[0]['mag_error'])
    logger.info("Photometry took: %.1f seconds", default_timer() - tic)

    return results_table


def do_plots(pm: PhotometryManager, mag_fig: plt.Figure, mag_ax: plt.Axes):
    # Plot the image:
    fig, ax = plt.subplots(1, 2, figsize=(20, 18))
    plot_image(pm.image.clean, ax=ax[0], scale='log', cbar='right', title='Image')
    plot_image(pm.image.mask, ax=ax[1], scale='linear', cbar='right', title='Mask')
    fig.savefig(pm.directories.save_as('original.png'), bbox_inches='tight')
    plt.close(fig)

    # Plot background estimation:
    fig, ax = plt.subplots(1, 3, figsize=(20, 6))
    plot_image(pm.image.clean, ax=ax[0], scale='log', title='Original')
    plot_image(pm.bkg.background, ax=ax[1], scale='log', title='Background')
    plot_image(pm.image.subclean, ax=ax[2], scale='log', title='Background subtracted')
    fig.savefig(pm.directories.save_as('background.png'), bbox_inches='tight')
    plt.close(fig)

    # Create plot of target and reference star positions:
    fig, ax = plt.subplots(1, 1, figsize=(20, 18))
    plot_image(pm.image.subclean, ax=ax, scale='log', cbar='right', title=pm.target.name)
    ax.scatter(pm.references.table['pixel_column'], pm.references.table['pixel_row'], c='r', marker='o', alpha=0.6)
    ax.scatter(pm.cleaner.gaussian_xys[:, 0], pm.cleaner.gaussian_xys[:, 1], marker='s', alpha=0.6, edgecolors='green',
               facecolors='none')
    ax.scatter(pm.target.pixel_column, pm.target.pixel_row, marker='+', s=20, c='r')
    fig.savefig(pm.directories.save_as('positions.png'), bbox_inches='tight')
    plt.close(fig)

    # Plot EPSF
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 15))
    plot_image(pm.psf_builder.epsf.data, ax=ax1, cmap='viridis')
    for a, ax in ((0, ax3), (1, ax2)):
        profile = pm.psf_builder.epsf.data.sum(axis=a)
        ax.plot(profile, 'k.-')
        ax.axvline(profile.argmax())
        ax.set_xlim(-0.5, len(profile) - 0.5)
    ax4.axis('off')
    fig.savefig(pm.directories.save_as('epsf.png'), bbox_inches='tight')
    plt.close(fig)
    del ax, ax1, ax2, ax3, ax4

    # Create two plots of the difference image:
    if pm.diff_im_exists:
        fig, ax = plt.subplots(1, 1, squeeze=True, figsize=(20, 20))
        plot_image(pm.diffimage.clean, ax=ax, cbar='right', title=pm.target.name)
        ax.plot(pm.target.pixel_column, pm.target.pixel_row, marker='+', markersize=20, color='r')
        fig.savefig(pm.directories.save_as('diffimg.png'), bbox_inches='tight')
        ax.set_xlim(pm.target.pixel_column - 50, pm.target.pixel_column + 50)
        ax.set_ylim(pm.target.pixel_row - 50, pm.target.pixel_row + 50)
        fig.savefig(pm.directories.save_as('diffimg_zoom.png'), bbox_inches='tight')
        plt.close(fig)

    # Calibration (handled in magnitudes.py).
    mag_fig.savefig(pm.directories.save_as('calibration.png'), bbox_inches='tight')
    plt.close(mag_fig)
