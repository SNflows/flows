# -*- coding: utf-8 -*-
"""

.. codeauthor:: Rasmus Handberg <rasmush@phys.au.dk>
"""

import argparse
import logging
from flows import api, download_catalog

if __name__ == '__main__':
    # Parse command line arguments:
    parser = argparse.ArgumentParser(description='Run catalog.')
    parser.add_argument('-d', '--debug', help='Print debug messages.', action='store_true')
    parser.add_argument('-q', '--quiet', help='Only report warnings and errors.', action='store_true')
    parser.add_argument('-t', '--target', type=str, help='Target to print catalog for.', nargs='?', default=None)
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

    for target in api.get_catalog_missing():
        logger.info("Downloading catalog for target=%s...", target)
        download_catalog(target)

    if args.target is not None:
        cat = api.get_catalog(args.target)

        print("Target:")
        cat['target'].pprint_all()
        print("\nReferences:")
        cat['references'].pprint_all()
        print("\nAvoid:")
        cat['avoid'].pprint_all()
