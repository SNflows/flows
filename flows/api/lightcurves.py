#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Fetch current lightcurve from Flows API.

.. codeauthor:: Rasmus Handberg <rasmush@phys.au.dk>
"""

import requests
import os.path
import tempfile
from astropy.table import Table
from ..config import load_config


# --------------------------------------------------------------------------------------------------
def get_lightcurve(target):
    """
    Retrieve lightcurve from Flows server.

    Parameters:
        target (int): Target to download lightcurve for.

    Returns:
        :class:`astropy.table.Table`: Table containing lightcurve.

    TODO:
        - Enable caching of files.

    .. codeauthor:: Rasmus Handberg <rasmush@phys.au.dk>
    """

    # Get API token from config file:
    config = load_config()
    token = config.get('api', 'token', fallback=None)
    if token is None:
        raise RuntimeError("No API token has been defined")

    # Send query to the Flows API:
    params = {'target': target}
    r = requests.get('https://flows.phys.au.dk/api/lightcurve.php', params=params,
                     headers={'Authorization': 'Bearer ' + token})
    r.raise_for_status()

    # Create tempory directory and save the file into there,
    # then open the file as a Table:
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpfile = os.path.join(tmpdir, 'table.ecsv')
        with open(tmpfile, 'w') as fid:
            fid.write(r.text)

        tab = Table.read(tmpfile, format='ascii.ecsv')

    return tab
