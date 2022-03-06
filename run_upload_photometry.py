#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Upload photometry results to Flows server.

.. codeauthor:: Rasmus Handberg <rasmush@phys.au.dk>
"""

import argparse
import logging
from flows import api


# --------------------------------------------------------------------------------------------------
def main():
    # Parse command line arguments:
    parser = argparse.ArgumentParser(description='Upload photometry.')
    parser.add_argument('-d', '--debug', help='Print debug messages.', action='store_true')
    parser.add_argument('-q', '--quiet', help='Only report warnings and errors.', action='store_true')
    parser.add_argument('fileids', type=int, help='File IDs to be uploaded.', nargs='+')
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
    logger = logging.getLogger('flows')
    if not logger.hasHandlers():
        logger.addHandler(console)
    logger.setLevel(logging_level)

    # Loop through the fileids and upload the results:
    for fid in args.fileids:
        api.upload_photometry(fid)


# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    main()
