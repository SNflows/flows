#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Clean bad source extraction, find and correct WCS.

.. codeauthor:: Emir Karamehmetoglu <emir.k@phys.au.dk>
.. codeauthor:: Rasmus Handberg <rasmush@phys.au.dk>
"""
from typing import Dict, Optional, TypeVar, Tuple, Union
from dataclasses import dataclass
import warnings

import numpy as np
from numpy.typing import ArrayLike, NDArray
import astroalign as aa
from astropy.coordinates import SkyCoord
from astropy import wcs
from astropy.stats import sigma_clip, gaussian_fwhm_to_sigma
from astropy.modeling import models, fitting
from astropy.time import Time
from astropy.utils.exceptions import ErfaWarning
import astropy.units as u
from astropy.table import Table
from copy import deepcopy
from bottleneck import nanmedian, nansum, nanmean, replace
from scipy.spatial import KDTree
import pandas as pd  # TODO: Convert to pure numpy implementation
#    From https://github.com/PJ-Watson/sep-pjw:
#    "The original release of sep by Kyle Barbary no longer appears to be maintained.
#    [...] The aim of sep-pjw is to offer a version of sep that resolves [support] issues,
#    whilst maintaining compatibility as much as is feasibly possible. [...]
#    For existing workflows, the only necessary update will be to change the import to [sep_pjw]"
#import sep
import sep_pjw


from .image import FlowsImage
from .target import Target
from .utilities import create_logger
logger = create_logger()

RefTable = TypeVar('RefTable', Dict, ArrayLike, Table)


class MinStarError(RuntimeError):
    pass


@dataclass
class References:
    table: RefTable
    coords: Optional[SkyCoord] = None
    mask: Optional[np.ndarray] = None  # positive mask ie True where we want it.
    xy: Optional[RefTable] = None

    def replace_nans_pm(self) -> None:
        replace(self.table['pm_ra'], np.NaN, 0.)
        replace(self.table['pm_dec'], np.NaN, 0.)

    def make_sky_coords(self, reference_time: float = 2015.5) -> None:
        self.replace_nans_pm()
        self.coords = SkyCoord(ra=self.table['ra'], dec=self.table['decl'], pm_ra_cosdec=self.table['pm_ra'],
                        pm_dec=self.table['pm_dec'], unit='deg', frame='icrs',
                        obstime=Time(reference_time, format='decimalyear'))

    def propagate(self, obstime: Time) -> None:
        if self.coords is None:
            raise AttributeError("References.coords is not defined.")
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", ErfaWarning)
            self.coords = self.coords.apply_space_motion(new_obstime=obstime)

    def copy(self) -> 'References':
        return References(self.__dataclass_fields__)

    @property
    def masked(self, mask: Optional[np.ndarray] = None) -> "References":
        if self.mask is None:
            if mask is None:
                raise AttributeError("No mask defined.")
            self.mask = mask

        copy = self.copy()
        for name in self.__dataclass_fields__:
            if getattr(self, name) is not None:
                setattr(copy, name, getattr(self, name)[self.mask])
        return copy

    def get_xy(self, img_wcs: wcs.WCS) -> None:
        """get pixel coordinates of reference stars"""
        self.xy = img_wcs.all_world2pix(list(zip(self.coords.ra.deg, self.coords.dec.deg)), 0)

    def make_pixel_columns(self, column_name: str ='pixel_column', row_name: str = 'pixel_row') -> None:
        self.table[column_name], self.table[row_name] = list(map(np.array, zip(*self.xy)))

    def _prepend_row(self, row: dict) -> None:
        self.table.insert_row(0, row)

    def add_target(self, target: Target, starid: int = 0) -> None:
        self._prepend_row(target.output_dict(starid=starid))
        if target.pixel_row and target.pixel_column:
            self.xy = np.vstack(((target.pixel_column, target.pixel_row), self.xy))


def use_sep(image: FlowsImage, tries: int = 5, thresh: float = 5.):

    # Use sep to for soure extraction
    sep_background = sep.Background(image.image, mask=image.mask)
    try:
        objects = sep.extract(image.image - sep_background, thresh=thresh, err=sep_background.globalrms,
                              mask=image.mask, deblend_cont=0.1, minarea=9, clean_param=2.0)
    except KeyboardInterrupt as e:
        raise e
    except Exception as e:
        logger.warning("SEP failed, trying again...")
        if tries > 0:
            thresh += 3
            return use_sep(image, tries - 1, thresh * 2)
        else:
            raise e
    sep_references = References(table=Table(objects))
    sep_references.xy = sep_references.table[['x', 'y']]
    sep_references.make_pixel_columns()
    return sep_references


def force_reject_g2d(xarray: ArrayLike, yarray: ArrayLike, image: Union[NDArray, np.ma.MaskedArray],
                     rsq_min: float = 0.5, radius: float = 10, fwhm_guess: float = 6.0, fwhm_min: float = 3.5,
                     fwhm_max: float = 18.0) -> Tuple[np.ma.MaskedArray, ...]:
    """
    It takes a list of x and y coordinates, and a 2D image, and returns a list of x and y coordinates,
    a list of r-squared values, and a boolean mask

    :param xarray: x-coordinates of the stars
    :type xarray: ArrayLike
    :param yarray: y-coordinates of the stars
    :type yarray: ArrayLike
    :param image: the image to be processed
    :type image: Union[NDArray, np.ma.MaskedArray]
    :param rsq_min: The minimum r-squared value for a star to be considered good
    :type rsq_min: float
    :param radius: The radius of the box around the star to fit, defaults to 10
    :type radius: float (optional)
    :param fwhm_guess: The initial guess for the FWHM of the star
    :type fwhm_guess: float
    :param fwhm_min: The minimum FWHM allowed for a star
    :type fwhm_min: float
    :param fwhm_max: The maximum FWHM allowed for a star to be considered good
    :type fwhm_max: float
    :return: masked_fwhms, masked_xys, mask, masked_rsqs
    """
    # Set up 2D Gaussian model for fitting to reference stars:
    g2d = models.Gaussian2D(amplitude=1.0, x_mean=radius, y_mean=radius, x_stddev=fwhm_guess * gaussian_fwhm_to_sigma)
    g2d.amplitude.bounds = (0.1, 2.0)
    g2d.x_mean.bounds = (0.5 * radius, 1.5 * radius)
    g2d.y_mean.bounds = (0.5 * radius, 1.5 * radius)
    g2d.x_stddev.bounds = (fwhm_min * gaussian_fwhm_to_sigma, fwhm_max * gaussian_fwhm_to_sigma)
    g2d.y_stddev.tied = lambda model: model.x_stddev
    g2d.theta.fixed = True
    gfitter = fitting.LevMarLSQFitter()

    # Stars reject
    N = len(xarray)
    fwhms = np.full((N, 2), np.NaN)
    xys = np.full((N, 2), np.NaN)
    rsqs = np.full(N, np.NaN)
    for i, (x, y) in enumerate(zip(xarray, yarray)):
        x = int(np.round(x))
        y = int(np.round(y))
        xmin = max(x - radius, 0)
        xmax = min(x + radius + 1, image.shape[1])
        ymin = max(y - radius, 0)
        ymax = min(y + radius + 1, image.shape[0])

        curr_star = deepcopy(image[ymin:ymax, xmin:xmax])

        edge = np.zeros_like(curr_star, dtype='bool')
        edge[(0, -1), :] = True
        edge[:, (0, -1)] = True
        curr_star -= nanmedian(curr_star[edge])
        curr_star /= np.nanmax(curr_star)

        ypos, xpos = np.indices(curr_star.shape)
        nan_filter = np.ones_like(curr_star, dtype='bool')
        nan_filter = nan_filter & np.isfinite(curr_star)
        if len(curr_star[nan_filter]) < 3:  # Not enough pixels to fit
            logger.debug(f"Not enough pixels to fit star, curr_star[nan_filter]:{curr_star[nan_filter]}")
            rsqs[i] = np.NaN
            fwhms[i] = np.NaN
            continue

        gfit = gfitter(g2d, x=xpos[nan_filter], y=ypos[nan_filter], z=curr_star[nan_filter])

        # Center
        xys[i] = np.array([gfit.x_mean + x - radius, gfit.y_mean + y - radius], dtype='float64')

        # Calculate rsq
        sstot = nansum((curr_star - nanmean(curr_star)) ** 2)
        sserr = nansum(gfitter.fit_info['fvec'] ** 2)
        rsqs[i] = 0 if sstot == 0 else 1.0 - (sserr / sstot)

        # FWHM
        fwhms[i] = gfit.x_fwhm

    masked_xys = np.ma.masked_array(xys, ~np.isfinite(xys))
    masked_rsqs = np.ma.masked_array(rsqs, ~np.isfinite(rsqs))
    mask = (masked_rsqs >= rsq_min) & (masked_rsqs < 1.0)  # Reject Rsq < rsq_min
    # changed
    # masked_xys = masked_xys[mask] # Clean extracted array.
    # to
    masked_xys.mask[~mask] = True
    # don't know if it breaks anything, but it doesn't make sense if
    # len(masked_xys) != len(masked_rsqs) FIXME
    masked_fwhms = np.ma.masked_array(fwhms, ~np.isfinite(fwhms))

    return masked_fwhms, masked_xys, mask.data, masked_rsqs


# --------------------------------------------------------------------------------------------------
def clean_with_rsq_and_get_fwhm(masked_fwhms, masked_rsqs, references, min_fwhm_references=2, min_references=6,
                                rsq_min=0.15):
    """
    Clean references and obtain fwhm using RSQ values.

    Parameters:
        masked_fwhms (np.ma.maskedarray): array of fwhms
        masked_rsqs (np.ma.maskedarray): array of rsq values
        references (astropy.table.Table): table of reference stars
        min_fwhm_references: (Default 2) min stars to get a fwhm
        min_references: (Default 6) min stars to aim for when cutting by R2
        rsq_min: (Default 0.15) min rsq value

    .. codeauthor:: Emir Karamehmetoglu <emir.k@phys.au.dk>
    """
    min_references_now = min_references
    rsqvals = np.arange(rsq_min, 0.95, 0.15)[::-1]
    fwhm_found = False
    min_references_achieved = False
    fwhm = np.nan
    # Clean based on R^2 Value
    while not min_references_achieved:
        for rsqval in rsqvals:
            mask = (masked_rsqs >= rsqval) & (masked_rsqs < 1.0)
            nreferences = np.sum(np.isfinite(masked_fwhms[mask]))
            if nreferences >= min_fwhm_references:
                _fwhms_cut_ = np.nanmean(sigma_clip(masked_fwhms[mask], maxiters=100, sigma=2.0))
                if not fwhm_found:
                    fwhm = _fwhms_cut_
                    fwhm_found = True
            if nreferences >= min_references_now:
                references = references[mask]
                min_references_achieved = True
                break
        if min_references_achieved:
            break
        min_references_now = min_references_now - 2
        if (min_references_now < 2) and fwhm_found:
            break
        elif not fwhm_found:
            raise RuntimeError("Could not estimate FWHM")

    if np.isnan(fwhm):
        raise RuntimeError("Could not estimate FWHM")

    # if minimum references not found, then take what we can get with even a weaker cut.
    # TODO: Is this right, or should we grab rsq_min (or even weaker?)
    min_references_now = min_references - 2
    while not min_references_achieved:
        mask = (masked_rsqs >= rsq_min) & (masked_rsqs < 1.0)
        nreferences = np.sum(np.isfinite(masked_fwhms[mask]))
        if nreferences >= min_references_now:
            references = references[mask]
            min_references_achieved = True
        rsq_min = rsq_min - 0.07
        min_references_now = min_references_now - 1

    # Check len of references as this is a destructive cleaning.
    # if len(references) == 2: logger.info('2 reference stars remaining, check WCS and image quality')
    if len(references) < 2:
        raise MinStarError(f"{len(references)} References remaining; could not estimate fwhm.")
    return fwhm, references


# --------------------------------------------------------------------------------------------------
def mkposxy(posx, posy):
    '''Make 2D np array for astroalign'''
    img_posxy = np.array([[x, y] for x, y in zip(posx, posy)], dtype="float64")
    return img_posxy


# --------------------------------------------------------------------------------------------------
def try_transform(source, target, pixeltol=2, nnearest=5, max_stars=50):
    aa.NUM_NEAREST_NEIGHBORS = nnearest
    aa.PIXEL_TOL = pixeltol
    transform, (sourcestars, targetstars) = aa.find_transform(source, target, max_control_points=max_stars)
    return sourcestars, targetstars


# --------------------------------------------------------------------------------------------------
def try_astroalign(source, target, pixeltol=2, nnearest=5, max_stars_n=50):
    # Get indexes of matched stars
    success = False
    try:
        source_stars, target_stars = try_transform(source, target, pixeltol=pixeltol, nnearest=nnearest,
                                                   max_stars=max_stars_n)
        source_ind = np.argwhere(np.in1d(source, source_stars)[::2]).flatten()
        target_ind = np.argwhere(np.in1d(target, target_stars)[::2]).flatten()
        success = True
    except aa.MaxIterError:
        source_ind, target_ind = 'None', 'None'
    return source_ind, target_ind, success


# --------------------------------------------------------------------------------------------------
def min_to_max_astroalign(source, target, fwhm=5, fwhm_min=1, fwhm_max=4, knn_min=5, knn_max=20, max_stars=100,
                          min_matches=3):
    """Try to find matches using astroalign asterisms by stepping through some parameters."""
    # Set max_control_points par based on number of stars and max_stars.
    nstars = max(len(source), len(source))
    if max_stars >= nstars:
        max_stars_list = 'None'
    else:
        if max_stars > 60:
            max_stars_list = (max_stars, 50, 4, 3)
        else:
            max_stars_list = (max_stars, 6, 4, 3)

    # Create max_stars step-through list if not given
    if max_stars_list == 'None':
        if nstars > 6:
            max_stars_list = (nstars, 5, 3)
        elif nstars > 3:
            max_stars_list = (nstars, 3)

    pixeltols = np.linspace(int(fwhm * fwhm_min), int(fwhm * fwhm_max), 4, dtype=int)
    nearest_neighbors = np.linspace(knn_min, min(knn_max, nstars), 4, dtype=int)

    success = False
    for max_stars_n in max_stars_list:
        for pixeltol in pixeltols:
            for nnearest in nearest_neighbors:
                source_ind, target_ind, success = try_astroalign(source, target, pixeltol=pixeltol, nnearest=nnearest,
                                                                 max_stars_n=max_stars_n)
                if success:
                    if len(source_ind) >= min_matches:
                        return source_ind, target_ind, success
                    else:
                        success = False
    return 'None', 'None', success


# --------------------------------------------------------------------------------------------------
def kdtree(source, target, fwhm=5, fwhm_max=4, min_matches=3):
    '''Use KDTree to get nearest neighbor matches within fwhm_max*fwhm distance'''

    # Use KDTree to rapidly efficiently query nearest neighbors

    tt = KDTree(target)
    st = KDTree(source)
    matches_list = st.query_ball_tree(tt, r=fwhm * fwhm_max)

    # indx = []
    targets = []
    sources = []
    for j, (sstar, match) in enumerate(zip(source, matches_list)):
        if np.array(target[match]).size != 0:
            targets.append(match[0])
            sources.append(j)
    sources = np.array(sources, dtype=int)
    targets = np.array(targets, dtype=int)

    # Return indexes of matches
    return sources, targets, len(sources) >= min_matches


# --------------------------------------------------------------------------------------------------
def get_new_wcs(extracted_ind, extracted_stars, clean_references, ref_ind, obstime, rakey='ra_obs', deckey='decl_obs'):
    targets = (extracted_stars[extracted_ind][:, 0], extracted_stars[extracted_ind][:, 1])

    c = SkyCoord(ra=clean_references[rakey][ref_ind], dec=clean_references[deckey][ref_ind], frame='icrs',
                 obstime=obstime)
    return wcs.utils.fit_wcs_from_points(targets, c)


def make_rsq_mask(masked_rsqs: np.ma.MaskedArray) -> np.ndarray:
    # Switching to pandas for easier selection
    df = pd.DataFrame(masked_rsqs, columns=['rsq'])
    return df.sort_values('rsq', ascending=False).dropna().index.values


def get_clean_references(reference_table: RefTable, masked_rsqs: np.ma.MaskedArray, min_references_ideal: int = 6,
                         min_references_abs: int = 3, rsq_min: float = 0.15, rsq_ideal: float = 0.5,
                         keep_max: int = 100, rescue_bad: bool = True) -> Tuple[RefTable, np.ndarray]:
    # Greedy first try
    mask = (masked_rsqs >= rsq_ideal) & (masked_rsqs < 1.0)
    mask = ~mask.data | mask.mask  # masked out of range values OR non-finite values
    masked_rsqs.mask = mask
    rsq_mask_index = make_rsq_mask(masked_rsqs)[:keep_max]
    if len(rsq_mask_index) >= min_references_ideal:
        return reference_table[rsq_mask_index], rsq_mask_index

    # Desperate second try
    mask = (masked_rsqs >= rsq_min) & (masked_rsqs < 1.0)
    mask = ~mask.data | mask.mask
    masked_rsqs.mask = mask
    rsq_mask_index = make_rsq_mask(masked_rsqs)[:min_references_ideal]
    if len(rsq_mask_index) >= min_references_abs:
        return reference_table[rsq_mask_index], rsq_mask_index
    if not rescue_bad:
        raise MinStarError(f'Less than {min_references_abs} clean stars and rescue_bad = False')

    # Extremely desperate last ditch attempt i.e. "rescue bad"
    mask = (masked_rsqs >= 0.02) & (masked_rsqs < 1.0)
    mask = ~mask.data | mask.mask
    masked_rsqs.mask = mask
    rsq_mask_index = make_rsq_mask(masked_rsqs)[:min_references_ideal]
    if len(rsq_mask_index) < 2:
        raise MinStarError('Less than 2 clean stars.')
    return reference_table[rsq_mask_index], rsq_mask_index  # Return if len >= 2


class ReferenceCleaner:

    def __init__(self, image: FlowsImage, references: References, rsq_min: float = 0.3,
                 min_references_ideal: int = 6, min_references_abs: int = 3):
        self.image = image
        self.references = references
        self.rsq_min = rsq_min
        self.min_references_ideal = min_references_ideal
        self.min_references_abs = min_references_abs
        self.gaussian_xys: Optional[np.ndarray] = None  # gaussian pixel positions

    def _clean_extracted_stars(self, x: Optional[ArrayLike] = None,
                               y: Optional[ArrayLike] = None) -> Tuple[np.ma.MaskedArray, ...]:
        """
        Clean extracted stars.
        :return: Tuple of masked_fwhms, masked_ref_xys, rsq_mask, masked_rsqs
        """
        # use instrument_defaults for initial guess of FWHM
        radius = self.image.instrument_defaults.radius
        fwhm_guess = self.image.instrument_defaults.fwhm
        fwhm_min = self.image.instrument_defaults.fwhm_min
        fwhm_max = self.image.instrument_defaults.fwhm_max
        useimage = self.image.subclean if self.image.subclean is not None else self.image.clean

        # Clean the references
        x = x if x is not None else self.references.table['pixel_column']
        y = y if y is not None else self.references.table['pixel_row']
        return force_reject_g2d(x, y, useimage, radius=radius, fwhm_guess=fwhm_guess, rsq_min=self.rsq_min,
                                fwhm_max=fwhm_max, fwhm_min=fwhm_min)

    def set_gaussian_xys(self, masked_ref_xys: np.ma.MaskedArray, old_references: RefTable,
                         new_references: RefTable) -> None:
        xy = [tuple(masked_ref_xys[old_references['starid'] == ref['starid']].data[0]) for ref in new_references]
        self.gaussian_xys = np.array(xy)

    def clean_references(self, references: References = None) -> Tuple[References, float]:
        if references is None:
            references = self.references

        # Clean the references
        masked_fwhms, masked_ref_xys, rsq_mask, masked_rsqs = self._clean_extracted_stars(
            references.table['pixel_column'],
            references.table['pixel_row'])


        # Use R^2 to more robustly determine initial FWHM guess.
        # This cleaning is good when we have FEW references.
        fwhm, fwhm_clean_references = clean_with_rsq_and_get_fwhm(
            masked_fwhms, masked_rsqs, references.table, min_fwhm_references=2,
            min_references=self.min_references_abs, rsq_min=self.rsq_min)
        logger.info('Initial FWHM guess is %f pixels', fwhm)

        # Final clean of wcs corrected references
        logger.info("Number of references before final cleaning: %d", len(references.table))
        logger.debug('Masked R^2 values: %s', masked_rsqs[rsq_mask].data)
        # Get references cleaned and ordered by R^2:
        ordered_cleaned_references, order_index = get_clean_references(references.table, masked_rsqs, rsq_ideal=0.8)
        ordered_coords = None if references.coords is None else references.coords[order_index]
        ordered_xy = None if references.xy is None else references.xy[order_index]
        ordered_cleaned_references = References(table=ordered_cleaned_references, coords=ordered_coords, xy=ordered_xy)
        logger.info("Number of references after final cleaning: %d", len(ordered_cleaned_references.table))

        # Save Gaussian XY positions before returning
        self.set_gaussian_xys(masked_ref_xys, references.table, ordered_cleaned_references.table)
        return ordered_cleaned_references, fwhm

    def make_sep_clean_references(self) -> References:
        """
        Make a clean reference catalog using SExtractor.
        """
        image = self.image
        # Get the SExtractor references from the image
        sep_references = use_sep(image)

        # Clean extracted stars
        _, masked_sep_xy, sep_mask, masked_sep_rsqs = self._clean_extracted_stars(
            sep_references.table['x'], sep_references.table['y'])

        sep_references.mask = sep_mask
        return sep_references.masked

    def mask_edge_and_target(self, target_coords: SkyCoord, hsize: int = 10,
                         target_distance_lim: u.quantity.Quantity = 10 * u.arcsec) -> References:
        """
        Clean the references by removing references that are too close to the target.
        """
        image_shape = self.image.shape

        # Make mask
        mask = (target_coords.separation(self.references.coords) > target_distance_lim) & (
            self.references.table['pixel_column'] > hsize) & (self.references.table['pixel_column'] < (image_shape[1] - 1 - hsize)) & (
            self.references.table['pixel_row'] > hsize) & (self.references.table['pixel_row'] < (image_shape[0] - 1 - hsize))
        self.references.mask = mask

        # Make new clean references
        return self.references.masked
        # return References(table=self.references.masked, mask=mask,
        #                  coords=self.references.coords[mask],
        #                  xy=self.references.xy[mask])


@dataclass(frozen=True)
class InitGuess:
    clean_references: References
    target_row: int = 0
    diff_row: Optional[int] = None

    @property
    def init_guess_full(self) -> Table:
        return Table(self.clean_references.xy, names=['x_0', 'y_0'])

    @property
    def init_guess_target(self) -> Table:
        return Table(self.init_guess_full[self.target_row])

    @property
    def init_guess_diff(self) -> Table:
        if self.diff_row is None:
            raise ValueError('`diff_row` is None, I cannot calculate the initial guesses for the difference image row.')
        return Table(self.init_guess_full[self.diff_row])

    @property
    def init_guess_references(self) -> Table:
        ref_begin = max(self.diff_row, self.target_row) + 1 if self.diff_row is not None else self.target_row + 1
        return self.init_guess_full[ref_begin:]
