#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Query TNS for new targets and upload to candidate marshal.
https://wis-tns.weizmann.ac.il/
TNS bot apikey must exist in config

.. codeauthor:: Emir Karamehmetoglu <emir.k@phys.au.dk>
.. codeauthor:: Rasmus Handberg <rasmush@phys.au.dk>
"""

import argparse
import logging
from tqdm import tqdm
import re
from astropy.coordinates import SkyCoord
from astropy.time import Time
from datetime import datetime, timedelta, timezone
from flows import api, tns


# --------------------------------------------------------------------------------------------------
def main():
    # Parse command line arguments:
    parser = argparse.ArgumentParser(description='Query TNS and upload to Flows candidates.')
    parser.add_argument('-d', '--debug', action='store_true', help='Print debug messages.')
    parser.add_argument('-q', '--quiet', action='store_true', help='Only report warnings and errors.')
    parser.add_argument('--zmax', type=float, default=0.105, help='Maximum redshift.')
    parser.add_argument('--zmin', type=float, default=0.000000001, help='Minimum redshift.')
    parser.add_argument('-b', '--days_begin', type=int, default=30, help='Discovery day at least X days before today.')
    parser.add_argument('-e', '--days_end', type=int, default=3, help='Discovery day at most X days before today.')
    parser.add_argument('-o', '--objtype', type=str, default=[3, 104],
                        help='TNS objtype int given as comma separed string with no spaces')
    parser.add_argument('-m', '--limit_months', type=int, default=2, help='Integer number of months to limit TNS search (for speed). \
                        Should be greater than days_begin.')
    parser.add_argument('--autoupload', action='store_true',
                        help="Automatically upload targets to Flows website. Only do this, if you know what you are doing!")
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

    tqdm_settings = {'disable': None if logger.isEnabledFor(logging.INFO) else True}

    # Try to load TNS config - only used for early stopping
    try:
        tns._load_tns_config()
    except tns.TNSConfigError:
        parser.error("Error in TNS configuration.")
        return

    # Calculate current date and date range to search:
    date_now = datetime.now(timezone.utc).date()
    date_end = date_now - timedelta(days=args.days_end)
    date_begin = date_now - timedelta(days=args.days_begin)
    logger.info('Date begin = %s, date_end = %s', date_begin, date_end)

    # Query TNS for SN names
    logger.info('Querying TNS for all targets, this may take awhile')
    nms = tns.tns_getnames(months=args.limit_months,  # pre-limit TNS search to candidates reported in the last X months
                           date_begin=date_begin, date_end=date_end, zmin=args.zmin, zmax=args.zmax,
                           objtype=args.objtype# Relevant TNS SN Ia subtypes.
                           )
    logger.debug(nms)

    if not nms:
        logger.info("No targets were found.")
        return

    # Remove already existing names using flows api
    included_names = ['SN' + target['target_name'] for target in api.get_targets()]
    nms = list(set(nms) - set(included_names))
    logger.info('Target names obtained: %s', nms)

    # Regular Expression matching any string starting with "ztf"
    regex_ztf = re.compile('^ztf', flags=re.IGNORECASE)
    regex_sn = re.compile(r'^sn\s*', flags=re.IGNORECASE)

    # Query TNS for object info using API, then upload to FLOWS using API.
    num_uploaded = 0
    if args.autoupload:
        for name in tqdm(nms, **tqdm_settings):
            sn = regex_sn.sub('', name)
            logger.debug('querying TNS for: %s', sn)

            # make GET request to TNS via API
            reply = tns.tns_get_obj(sn)

            # Parse output
            if reply:
                logger.debug('GET query successful')

                # Extract object info
                coord = SkyCoord(ra=reply['radeg'], dec=reply['decdeg'], unit='deg', frame='icrs')
                discovery_date = Time(reply['discoverydate'], format='iso', scale='utc')
                ztf = list(filter(regex_ztf.match, reply['internal_names']))
                ztf = None if not ztf else ztf[0]
                if 'object_type' in reply and 'name' in reply['object_type']:
                    sntype = regex_sn.sub('', reply['object_type']['name'])
                else:
                    sntype = None

                # Try to upload to FLOWS
                newtargetid = api.add_target(reply['objname'], coord, redshift=reply['redshift'],
                                             discovery_date=discovery_date, discovery_mag=reply['discoverymag'],
                                             host_galaxy=reply['hostname'], ztf=ztf, sntype=sntype, status='candidate',
                                             project='flows')
                logger.debug('Uploaded to FLOWS with targetid=%d', newtargetid)
                num_uploaded += 1

    logger.info("%d targets uploaded to Flows.", num_uploaded)


# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    main()
