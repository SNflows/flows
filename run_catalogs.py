"""
Runner to add target catalog, if catalog info is missing. Print output.
"""
import argparse
import logging
from tendrils import api
from flows import download_catalog


def parse():
    """
     # Parse command line arguments:
    """
    parser = argparse.ArgumentParser(description='Run catalog.')
    parser.add_argument('-d', '--debug', help='Print debug messages.', action='store_true')
    parser.add_argument('-q', '--quiet', help='Only report warnings and errors.', action='store_true')
    parser.add_argument('-c', '--commit', help='Commit downloaded catalog(s) to flows database.', action='store_true')
    parser.add_argument('-p', '--print', help='Print single catalog which already exists in flows database.', action='store_true')
    parser.add_argument('-t', '--target', type=str, help='Optionally specify just one target for downloading/committing/printing.', nargs='?', default=None)
    return parser.parse_args()

def set_logging_level(args):
    # Set logging level:
    logging_level = logging.INFO
    if args.quiet:
        logging_level = logging.WARNING
    elif args.debug:
        logging_level = logging.DEBUG
    return logging_level

def main():
    # Parse command line arguments:
    args = parse()

    # Setup logging:
    logging_level = set_logging_level(args)
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    console = logging.StreamHandler()
    console.setFormatter(formatter)
    logger = logging.getLogger('flows')
    if not logger.hasHandlers():
        logger.addHandler(console)
    logger.setLevel(logging_level)

    targets = api.get_catalog_missing()
    if args.target is not None:
        if args.print is None:
            # We want to download and possibly commit the target catalog
            if int(args.target) in targets:
                logger.info("Downloading catalog for target=%s (committing to db=%s)...", args.target, args.commit)
                download_catalog(args.target, update_existing=args.commit)  # @TODO: refactor to Tendrils
            else:
                logger.warning("Cannot find target=%s in list generated by api.get_catalog_missing()", args.target)
        else:
            # Download target catalog from db for printing
            cat = api.get_catalog(args.target)
            print(f"Target:{cat['target'].pprint_all()} "
                  f"\nReferences: {cat['references'].pprint_all()} "
                  f"\nAvoid:cat['avoid'].pprint_all()")
    else:
        # Download and possibly commit all missing target catalogs
        logger.info("%d catalogs missing", len(targets))
        for target in targets:
            logger.info("Downloading catalog for target=%s...", target)
            download_catalog(target, update_existing=args.commit)  # @TODO: refactor to Tendrils

if __name__ == '__main__':
    main()
