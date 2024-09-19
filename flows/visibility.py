#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Target visibility plotting.
@TODO: Move to flows-tools.

.. codeauthor:: Rasmus Handberg <rasmush@phys.au.dk>
"""

import logging
import os.path
import warnings
from datetime import datetime
from typing import Optional

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from astropy.coordinates import AltAz, SkyCoord, get_body, get_sun
from astropy.time import Time
from astropy.visualization import quantity_support
from matplotlib import pyplot
from matplotlib.dates import DateFormatter
from tendrils import api

from .target import Target


# --------------------------------------------------------------------------------------------------
def visibility(target: Target, siteid: Optional[int] = None, date=None, output=None, overwrite=True):
    """
    Create visibility plot.

    Parameters:
        target (str or int):
        siteid (int): Identifier of site.
        date (datetime or str, optional): Date for which to create visibility plot.
            Default it to use the current date.
        output (str, optional): Path to file or directory where to place the plot.
            If not given, the plot will be created in memory, and can be shown on screen.
        overwrite (bool, optional): Should existing file specified in ``output`` be overwritten?
            Default is to overwrite an existing file.

    .. codeauthor:: Rasmus Handberg <rasmush@phys.au.dk>
    """

    logger = logging.getLogger(__name__)
    warnings.warn(DeprecationWarning("This module is moved to SNFLOWS/flows-tools."))

    if date is None:
        date = datetime.utcnow()
    elif isinstance(date, str):
        date = datetime.strptime(date, '%Y-%m-%d')

    # Coordinates of object:
    obj = SkyCoord(ra=target.ra, dec=target.dec, unit='deg', frame='icrs')

    if siteid is None:
        sites = api.get_all_sites()
    else:
        sites = [api.get_site(siteid)]

    plotpaths = []
    for site in sites:
        # If we are saving plot to file, determine the path to save to
        # and check that it doesn't already exist:
        if output:
            if os.path.isdir(output):
                plotpath = os.path.join(output, "visibility_%s_%s_site%02d.png" % (
                    target.name, date.strftime('%Y%m%d'), site['siteid']))
            else:
                plotpath = output
            logger.debug("Will save visibility plot to '%s'", plotpath)

            # If we are not overwriting and
            if not overwrite and os.path.exists(plotpath):
                logger.info("File already exists: %s", plotpath)
                continue
        logger.debug(f"{site}, type: {type(site)}")
        # Observatory:
        observatory = site.get('EarthLocation', None)
        if observatory is None:
            observatory = site.get('earth_location')
        utcoffset = (site['longitude'] * u.deg / (360 * u.deg)) * 24 * u.hour

        # Create timestamps to calculate for:
        midnight = Time(date.strftime('%Y-%m-%d') + ' 00:00:00', scale='utc') - utcoffset
        delta_midnight = np.linspace(-12, 12, 1000) * u.hour
        times = midnight + delta_midnight

        # AltAz frame:
        AltAzFrame = AltAz(obstime=times, location=observatory)

        # Object:
        altaz_obj = obj.transform_to(AltAzFrame)

        # The Sun and Moon:
        altaz_sun = get_sun(times).transform_to(AltAzFrame)
        altaz_moon = get_body("moon", times).transform_to(AltAzFrame)

        sundown_astro = (altaz_sun.alt < -6 * u.deg)
        if np.any(sundown_astro):
            min_time = np.min(times[sundown_astro]) - 2 * u.hour
            max_time = np.max(times[sundown_astro]) + 2 * u.hour
        else:
            min_time = times[0]
            max_time = times[-1]

        quantity_support()
        fig, ax = plt.subplots(1, 1, figsize=(15, 9), squeeze=True)
        plt.grid(ls=':', lw=0.5)
        ax.plot(times.datetime, altaz_sun.alt, color='y', label='Sun')
        ax.plot(times.datetime, altaz_moon.alt, color=[0.75] * 3, ls='--', label='Moon')
        objsc = ax.scatter(times.datetime, altaz_obj.alt.deg, c=altaz_obj.az.deg, label=target.name, lw=0, s=8,
                           cmap='twilight')
        ax.fill_between(times.datetime, 0 * u.deg, 90 * u.deg, altaz_sun.alt < -0 * u.deg, color='0.5',
                        zorder=0)  # , label='Night'
        ax.fill_between(times.datetime, 0 * u.deg, 90 * u.deg, altaz_sun.alt < -18 * u.deg, color='k',
                        zorder=0)  # , label='Astronomical Night'

        plt.colorbar(objsc, ax=ax, pad=0.01).set_label('Azimuth [deg]')
        ax.legend(loc='upper left')
        ax.minorticks_on()
        ax.set_xlim(min_time.datetime, max_time.datetime)
        ax.set_ylim(0 * u.deg, 90 * u.deg)
        ax.set_title("%s - %s - %s" % (target.name, date.strftime('%Y-%m-%d'), site['sitename']),
                     fontsize=14)
        plt.xlabel('Time [UTC]', fontsize=14)
        plt.ylabel('Altitude [deg]', fontsize=16)
        fig.autofmt_xdate()

        formatter = DateFormatter('%d/%m %H:%M')
        ax.xaxis.set_major_formatter(formatter)

        if output:
            fig.savefig(plotpath, bbox_inches='tight', transparent=True)
            plt.close(fig)
            plotpaths.append(plotpath)

    if output:
        return plotpaths

    if pyplot.isinteractive():
        plt.show()
    return ax
