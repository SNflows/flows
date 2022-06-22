#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# flake8: noqa

import numpy as np
import astropy.wcs
from .coordinatematch import CoordinateMatch
from .wcs import WCS2
from flows.utilities import create_logger
from astropy.coordinates import SkyCoord
from flows.image import FlowsImage
from flows import reference_cleaning as refclean
from flows.target import Target

logger = create_logger()


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
        image.wcs = astropy.wcs.utils.fit_wcs_from_points(np.array(list(zip(*cm.xy[i_xy]))),
                                                          SkyCoord(*map(list, zip(*cm.rd[i_rd])), unit='deg'))
        del i_xy, i_rd

    logger.debug(f'Used WCS: {WCS2.from_astropy_wcs(image.wcs)}')
    return image
