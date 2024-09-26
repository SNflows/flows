#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

.. codeauthor:: Simon Holmbo
.. codeauthor:: Rasmus Handberg <rasmush@phys.au.dk>
"""

import os.path
import sys
import tempfile
import shutil
import subprocess
import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from astropy.io import fits
from astropy.coordinates import SkyCoord
import healpy
if os.path.abspath(os.path.join(os.path.dirname(__file__), 'flows')) not in sys.path:
    sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), 'flows')))
from tendrils import api
from tendrils.utils import urls_from_config
from flows.plots import plots_noninteractive
import urllib.request

#--------------------------------------------------------------------------------------------------
def optipng(fpath):
    subprocess.run(['optipng', '-preserve', '-quiet', fpath])

#--------------------------------------------------------------------------------------------------
def healpix_to_cartesian(data):
    fig = plt.figure(0)
    healpy.cartview(data, fig=0, nest=True)
    data = fig.axes[0].images[0].get_array()
    fig.clf()
    plt.close(fig)
    return data

#--------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    # Shorthands used later on:
    __dir__ = os.path.dirname(os.path.abspath(__file__))

    # Switch on non-interactive backend:
    plots_noninteractive()

    # download milkyway data if it does not exist locally
    filename = os.path.join(__dir__, "flows_mwmap_milkyway.fits")
    if not os.path.exists(filename):
        mwfile = urllib.request.urlretrieve(url='https://anon.erda.au.dk/share_redirect/fnsq9iKMmp', filename=filename)

    # Load milky-way data:
    print("loading mw data:", filename)
    with fits.open(filename, mode='readonly') as hdul:
        data = hdul[1].data["temperature"]

    data = healpix_to_cartesian(data)
    LON, LAT = np.meshgrid(np.linspace(180, -180, data.shape[1]) % 360, np.linspace(-90, 90, data.shape[0]))
    RADEC = SkyCoord(LON, LAT, frame="galactic", unit="degree")
    RA, DEC = RADEC.icrs.ra.degree, RADEC.icrs.dec.degree
    ra, dec = np.meshgrid(np.linspace(180, -180, data.shape[1]) % 360, np.linspace(-90, 90, data.shape[0]))
    data = interpolate.griddata(list(zip(RA.ravel(), DEC.ravel())), data.ravel(), list(zip(ra.ravel(), dec.ravel()))).reshape(data.shape)

    # FLOWS
    urls = urls_from_config()
    print("loading flows data:", urls.targets_url)
    dtype = {"names": ("name", "ra", "dec", "redshift"), "formats": ("U10", 'f', 'f', 'f')}
    targets = np.array([(obj["target_name"], obj["ra"], obj["decl"], obj["redshift"]) for obj in api.get_targets() if obj["project"] == "flows" and obj["target_status"] == "target"], dtype=dtype)

    print("plotting data...", np.size(targets, axis=0), "targets")
    fig = plt.figure(figsize=(8, 4))
    ax = plt.subplot(projection="mollweide")
    X, Y = np.meshgrid(np.linspace(np.pi, -np.pi, data.shape[1]), np.linspace(-np.pi/2, np.pi/2, data.shape[0]))
    vmin, vmax = np.nanquantile(np.log(data), (.01, .99))
    ax.pcolormesh(X, Y, np.log(data), cmap="gray_r", vmin=vmin, vmax=vmax, shading="auto", rasterized=True)

    ra = (np.deg2rad(targets["ra"]) + np.pi) % (2*np.pi) - np.pi
    dec = np.deg2rad(targets["dec"])

    vmin, vmax = 0.00, 0.10 # np.nanquantile(targets["redshift"], (.01, .99))
    points = ax.scatter(ra, dec, s=10, c=targets["redshift"], cmap="Reds", vmin=vmin, vmax=vmax, zorder=2)

    for x in (-150, -120, -90, -60, -30, 0, 30, 60, 90, 120, 150):
        ax.text(np.deg2rad(x), 0, fr"{x}$\degree$", color="black", ha="center", va="center", zorder=3)
        ax.add_patch(Rectangle((np.deg2rad(x)-.2, -.05), .4, .1325, color="white", linewidth=0, alpha=.65, zorder=2))
    ax.set_xticklabels([])

    ax.grid(color="black", ls="dotted", lw=.5, zorder=1)

    N = 1000

    cbar = ax.inset_axes([0.5, 0.0, 0.6, 1.0], zorder=-1)
    X, Y = np.meshgrid(*map(np.arange, (N, N)))
    Z = -np.arctan2(Y - N/2, X - N/2)
    Z[Z <= -np.pi/4] = np.nan
    Z[Z >= +np.pi/4] = np.nan
    Z[(X - N/2)**2 + (Y - N/2)**2 <= (N/2.9)**2] = np.nan
    Z[(X - N/2)**2 + (Y - N/2)**2 >= (N/2.5)**2] = np.nan
    cbar.imshow(Z, cmap="Reds", aspect="auto", interpolation="nearest")

    xy = [N/2 + N/2.5 * np.array([np.cos(t), np.sin(t)]) for t in np.linspace(-np.pi/4, np.pi/4, N//5)]
    xy += [N/2 + N/2.9 * np.array([np.cos(t), np.sin(t)]) for t in np.linspace(+np.pi/4, -np.pi/4, N//5)] + [xy[0]]
    cbar.plot(*zip(*xy), marker=None, linewidth=1, color="black")

    for i, t in enumerate(np.linspace(-np.pi/4, +np.pi/4, 9+10+2)[1:-1]):
        s = 0.1 if not (i + 1) % 2 else 0.05
        xy  = [N/2 + N/2.5 * np.array([np.cos(t), np.sin(t)])]
        xy += [N/2 + N/(2.5+s) * np.array([np.cos(t), np.sin(t)])]
        cbar.plot(*zip(*xy), marker=None, linewidth=1, color="black")

    cbar.text(N/2 + N/2.7 * np.cos(+np.pi/3.25), N/2 + N/2.7 * np.sin(+np.pi/3.25), f"z = {vmin:.2f}", va="top")
    cbar.text(N/2 + N/2.7 * np.cos(-np.pi/3.25), N/2 + N/2.7 * np.sin(-np.pi/3.25), f"z = {vmax:.2f}", va="bottom")
    cbar.axis("off")

    # Work in temp directory:
    with tempfile.TemporaryDirectory() as tmpdir:
        # Save figure:
        fpath = os.path.join(tmpdir, 'mwmap.png')
        fig.savefig(fpath, bbox_inches='tight', transparent=True)

        # Optimize the generated image:
        optipng(fpath)

        # Overwrite existing file:
        shutil.move(fpath, os.path.join(__dir__, "mwmap.png"))
