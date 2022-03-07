#!/usr/bin/env python3
"""
Download ZTF photometry from ALERCE API.
https://alerceapi.readthedocs.io/
"""

import argparse
import logging
import os
import numpy as np
import matplotlib.pyplot as plt
from tendrils import api
from tendrils.utils import load_config, ztf


def main():
    # Parse command line arguments:
    parser = argparse.ArgumentParser(description='Download ZTF photometry.')
    parser.add_argument('-d', '--debug', help='Print debug messages.', action='store_true')
    parser.add_argument('-q', '--quiet', help='Only report warnings and errors.', action='store_true')
    parser.add_argument('-t', '--target', type=str, default=None, help='Target to download ZTF photometry for.')
    parser.add_argument('-o', '--output', type=str, default=None, help='Directory to save output to.')
    args = parser.parse_args()

    # Set logging level:
    logging_level = logging.INFO
    if args.quiet:
        logging_level = logging.WARNING
    elif args.debug:
        logging_level = logging.DEBUG

    # Setup logging:
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    console = logging.StreamHandler()
    console.setFormatter(formatter)
    logger = logging.getLogger(__name__)
    if not logger.hasHandlers():
        logger.addHandler(console)
    logger.setLevel(logging_level)

    if args.output is None:
        config = load_config()
        output_dir = config.get('ztf', 'output_photometry', fallback='.')
    else:
        output_dir = args.output
    logger.info("Saving output to '%s'", output_dir)

    # Check that output directory exists:
    if not os.path.isdir(output_dir):
        parser.error(f"Output directory does not exist: '{output_dir}'")  # noqa: G004

    # Use API to get list of targets to process:
    if args.target is None:
        targets = api.get_targets()
    else:
        targets = [api.get_target(args.target)]

    # Colors used for the different filters in plots:
    # I know purple is in the wrong end of the scale, but not much I can do
    colors = {'gp': 'tab:green', 'rp': 'tab:red', 'ip': 'tab:purple'}

    # Loop through targets:
    for tgt in targets:
        logger.debug("Target: %s", tgt)
        target_name = tgt['target_name']

        # Paths to the files to be updated:
        ztf_lightcurve_path = os.path.join(output_dir, f'{target_name:s}-ztf.ecsv')
        ztf_plot_path = os.path.join(output_dir, f'{target_name:s}-ztf.png')

        # If there is no ZTF id, there is no need to try:
        # If an old file exists then delete it.
        if tgt['ztf_id'] is None:
            if os.path.isfile(ztf_lightcurve_path):
                os.remove(ztf_lightcurve_path)
            if os.path.isfile(ztf_plot_path):
                os.remove(ztf_plot_path)
            continue

        # Download ZTF photometry as Astropy Table:
        tab = ztf.download_ztf_photometry(tgt['targetid'])
        logger.debug("ZTF Photometry:\n%s", tab)
        if tab is None or len(tab) == 0:
            if os.path.isfile(ztf_lightcurve_path):
                os.remove(ztf_lightcurve_path)
            if os.path.isfile(ztf_plot_path):
                os.remove(ztf_plot_path)
            continue

        # Write table to file:
        tab.write(ztf_lightcurve_path, format='ascii.ecsv', delimiter=',')

        # Find time of maxmimum and 14 days from that:
        indx_min = np.argmin(tab['mag'])
        maximum_mjd = tab['time'][indx_min]
        fortnight_mjd = maximum_mjd + 14

        # Get LC data out and save as CSV files
        fig, ax = plt.subplots()
        ax.axvline(maximum_mjd, ls='--', c='k', lw=0.5, label='Maximum')
        ax.axvline(fortnight_mjd, ls='--', c='0.5', lw=0.5, label='+14 days')
        for fid in np.unique(tab['photfilter']):
            col = colors[fid]
            band = tab[tab['photfilter'] == fid]
            ax.errorbar(band['time'], band['mag'], band['mag_err'], color=col, ls='-', lw=0.5, marker='.', label=fid)

        ax.invert_yaxis()
        ax.set_title(target_name)
        ax.set_xlabel('Time (MJD)')
        ax.set_ylabel('Magnitude')
        ax.legend()
        fig.savefig(ztf_plot_path, format='png', bbox_inches='tight')
        plt.close(fig)


# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    main()
