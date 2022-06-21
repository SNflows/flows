"""
Run Flows photometry. Allows multithreaded operations to be run
"""

import argparse
import logging
import os
import sys
import shutil
import functools
import multiprocessing
import tqdm
from tendrils import api, utils
from flows import photometry


def process_fileid(fid, output_folder_root=None, attempt_imagematch=True, autoupload=False, keep_diff_fixed=False,
                   cm_timeout=None, noplots=False):
    logger = logging.getLogger('flows')
    logging.captureWarnings(True)
    logger_warn = logging.getLogger('py.warnings')
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s', "%Y-%m-%d %H:%M:%S")

    datafile = api.get_datafile(fid)
    target_name = datafile['target_name']

    # Folder to save output:
    output_folder = os.path.join(output_folder_root, target_name, f'{fid:05d}')

    photfile = None
    _filehandler = None
    try:
        # Set the status to indicate that we have started processing:
        if autoupload:
            api.set_photometry_status(fid, 'running')

        # Create the output directory if it doesn't exist:
        os.makedirs(output_folder, exist_ok=True)

        # Also write any logging output to the
        _filehandler = logging.FileHandler(os.path.join(output_folder, 'photometry.log'), mode='w')
        _filehandler.setFormatter(formatter)
        _filehandler.setLevel(logging.INFO)
        logger.addHandler(_filehandler)
        logger_warn.addHandler(_filehandler)

        photfile = photometry(fileid=fid, cm_timeout=cm_timeout, make_plots=not noplots)
        # photfile = photometry(fileid=fid, output_folder=output_folder, attempt_imagematch=attempt_imagematch,
        #                       keep_diff_fixed=keep_diff_fixed, cm_timeout=cm_timeout)

    except (SystemExit, KeyboardInterrupt):
        logger.error("Aborted by user or system.")
        if os.path.exists(output_folder):
            shutil.rmtree(output_folder, ignore_errors=True)
        photfile = None
        if autoupload:
            api.set_photometry_status(fid, 'abort')

    except:  # noqa: E722, pragma: no cover
        logger.exception("Photometry failed")
        photfile = None
        if autoupload:
            api.set_photometry_status(fid, 'error')

    if _filehandler is not None:
        logger.removeHandler(_filehandler)
        logger_warn.removeHandler(_filehandler)

    if photfile is not None:
        if autoupload:
            api.upload_photometry(fid, delete_completed=True)
            api.set_photometry_status(fid, 'ingest')

    return photfile


def main():
    # Parse command line arguments:
    parser = argparse.ArgumentParser(description='Run photometry pipeline.')
    parser.add_argument('-d', '--debug', help='Print debug messages.', action='store_true')
    parser.add_argument('-q', '--quiet', help='Only report warnings and errors.', action='store_true')
    parser.add_argument('-o', '--overwrite', help='Overwrite existing results.', action='store_true')

    group = parser.add_argument_group('Selecting which files to process')
    group.add_argument('--fileid', type=int, default=None, action='append',
                       help="Process this file ID. Overrides all other filters.")
    group.add_argument('--targetid', type=int, default=None, action='append',
                       help="Only process files from this target.")
    group.add_argument('--filter', type=str, default=None, choices=['missing', 'all', 'error'])
    group.add_argument('--minversion', type=str, default=None,
                       help="Include files not previously processed with at least this version.")

    group = parser.add_argument_group('Processing settings')
    group.add_argument('--threads', type=int, default=1, help="Number of parallel threads to use.")
    group.add_argument('--noplots', action='store_true', help="Disable plots")
    group.add_argument('--no-imagematch', action='store_true', help="Disable ImageMatch.")
    group.add_argument('--autoupload', action='store_true',
                       help="Automatically upload completed photometry to Flows website. Only do this, if you know what you are doing!")
    group.add_argument('--fixposdiff', action='store_true',
                       help="Fix SN position during PSF photometry of difference image. Useful when difference image is noisy.")
    group.add_argument('--wcstimeout', type=int, default=None, help="Timeout in Seconds for WCS.")
    args = parser.parse_args()

    # Ensure that all input has been given:
    if not args.fileid and not args.targetid and args.filter is None:
        parser.error("Please select either a specific FILEID .")

    # Set logging level:
    logging_level = logging.INFO
    if args.quiet:
        logging_level = logging.WARNING
    elif args.debug:
        logging_level = logging.DEBUG

    # Number of threads to use:
    threads = args.threads
    if threads <= 0:
        threads = multiprocessing.cpu_count()

    # Setup logging:
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s', "%Y-%m-%d %H:%M:%S")
    console = logging.StreamHandler(sys.stdout)
    console.setFormatter(formatter)
    logger = logging.getLogger('flows')
    if not logger.hasHandlers():
        logger.addHandler(console)
    logger.propagate = False
    logger.setLevel(logging_level)

    # If we have started a new processing, perform a cleanup of the
    # photometry status indicator. This will change all processes
    # still marked as "running" to "abort" if they have been running
    # for more than a day:
    if args.autoupload:
        api.cleanup_photometry_status()

    if args.fileid is not None:
        # Run the specified fileids:
        fileids = args.fileid
    else:
        # Ask the API for a list of fileids which are yet to be processed:
        if args.targetid is not None:
            fileids = []
            for targid in args.targetid:
                fileids += api.get_datafiles(targetid=targid, filt=args.filter, minversion=args.minversion)
        else:
            fileids = api.get_datafiles(filt=args.filter, minversion=args.minversion)

    # Remove duplicates from fileids to be processed:
    fileids = list(set(fileids))

    # Ask the config where we should store the output:
    config = utils.load_config()
    output_folder_root = config.get('photometry', 'output', fallback='.')

    # Create function wrapper:
    process_fileid_wrapper = functools.partial(process_fileid, output_folder_root=output_folder_root,
                                               attempt_imagematch=not args.no_imagematch, autoupload=args.autoupload,
                                               keep_diff_fixed=args.fixposdiff, cm_timeout=args.wcstimeout,
                                               noplots=args.noplots)

    if threads > 1:
        # Disable printing info messages from the parent function.
        # It is going to be all jumbled up anyway.
        if not args.debug:
            logger.setLevel(logging.WARNING)  # Don't silence if debugging.

        # There is more than one area to process, so let's start
        # a process pool and process them in parallel:
        with multiprocessing.Pool(threads) as pool:
            for result in tqdm.tqdm(pool.imap_unordered(process_fileid_wrapper, fileids), total=len(fileids)):
                # We can do something with the partial results here (like calculate color terms!).
                pass

    else:
        # Only single thread so simply run it directly:
        for fid in fileids:
            logger.info("=" * 72)
            logger.info(fid)
            logger.info("=" * 72)
            process_fileid_wrapper(fid)


if __name__ == '__main__':
    main()
