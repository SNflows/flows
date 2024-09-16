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
    parser.add_argument('-t', '--target', type=str, help='Target to print catalog for.', nargs='?', default=None)
    parser.add_argument('-c', '--commit', help='Commit downloaded catalogs to flows database.', action='store_true')
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

    # Get missing
    for target in api.get_catalog_missing():
        logger.info("Downloading catalog for target=%s...", target)
        download_catalog(target, update_existing=args.commit)  # @TODO: refactor to Tendrils

    # download target catalog for printing
    if args.target is not None:
        cat = api.get_catalog(args.target)

        print(f"Target:{cat['target'].pprint_all()} "
              f"\nReferences: {cat['references'].pprint_all()} "
              f"\nAvoid:cat['avoid'].pprint_all()")


if __name__ == '__main__':
    main()
