#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

.. codeauthor:: Rasmus Handberg <rasmush@phys.au.dk>
"""

import logging
import requests
from ..config import load_config


# --------------------------------------------------------------------------------------------------
def set_photometry_status(fileid, status):
    """
    Set photometry status.

    Parameters:
        fileid (int):
        status (str): Choises are 'running', 'error' or 'done'.

    .. codeauthor:: Rasmus Handberg <rasmush@phys.au.dk>
    """
    # Validate the input:
    logger = logging.getLogger(__name__)
    if status not in ('running', 'error', 'abort', 'ingest', 'done'):
        raise ValueError('Invalid status')

    # Get API token from config file:
    config = load_config()
    i_am_pipeline = config.getboolean('api', 'pipeline', fallback=False)
    if not i_am_pipeline:
        logger.debug("Not setting status since user is not pipeline")
        return False

    # Get API token from config file:
    token = config.get('api', 'token', fallback=None)
    if token is None:
        raise RuntimeError("No API token has been defined")

    # Send HTTP request to FLOWS server:
    r = requests.get('https://flows.phys.au.dk/api/set_photometry_status.php',
                     params={'fileid': fileid, 'status': status}, headers={'Authorization': 'Bearer ' + token})
    r.raise_for_status()
    res = r.text.strip()

    if res != 'OK':
        raise RuntimeError(res)

    return True


# --------------------------------------------------------------------------------------------------
def cleanup_photometry_status():
    """
    Perform a cleanup of the photometry status indicator.

    This will change all processes still marked as "running"
    to "abort" if they have been running for more than a day.

    .. codeauthor:: Rasmus Handberg <rasmush@phys.au.dk>
    """
    # Validate the input:
    logger = logging.getLogger(__name__)

    # Get API token from config file:
    config = load_config()
    i_am_pipeline = config.getboolean('api', 'pipeline', fallback=False)
    if not i_am_pipeline:
        logger.debug("Not setting status since user is not pipeline")
        return False

    # Get API token from config file:
    token = config.get('api', 'token', fallback=None)
    if token is None:
        raise RuntimeError("No API token has been defined")

    # Send HTTP request to FLOWS server:
    r = requests.get('https://flows.phys.au.dk/api/cleanup_photometry_status.php',
                     headers={'Authorization': 'Bearer ' + token})
    r.raise_for_status()
    res = r.text.strip()

    if res != 'OK':
        raise RuntimeError(res)

    return True
