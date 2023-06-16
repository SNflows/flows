"""
Run Flows photometry. Allows multithreaded operations to be run
"""

import argparse
import functools
import multiprocessing
import tqdm
from tendrils import api
from flows import photometry, fileio, result_model
from flows.utilities import create_logger, parse_log_level, create_warning_logger, remove_file_handlers


def process_fileid(fid, autoupload: bool = False, cm_timeout=None, no_plots: bool = False,
                   rescale: bool = True, rescale_dynamic: bool = True) -> result_model.ResultsTable:
    # Create the output directory if it doesn't exist:
    datafile = api.get_datafile(fid)
    directories = fileio.Directories.from_fid(fid, datafile=datafile)
    logger = create_logger(str(fid), log_file=directories.log_path)
    logger_warn = create_warning_logger(log_file=directories.log_path)

    try:
        # Set the status to indicate that we have started processing:
        if autoupload:
            api.set_photometry_status(fid, 'running')

        table = photometry(fileid=fid, cm_timeout=cm_timeout, make_plots=not no_plots,
                           directories=directories, datafile=datafile, rescale=rescale, 
                           rescale_dynamic=rescale_dynamic)

    except (SystemExit, KeyboardInterrupt):
        logger.error("Aborted by user or system.")
        fileio.del_dir(directories.output_folder)
        table = None
        if autoupload:
            api.set_photometry_status(fid, 'abort')

    except Exception as e:  # noqa: E722, pragma: no cover
        logger.exception(f"Photometry failed:{str(e)}")
        table = None
        if autoupload:
            api.set_photometry_status(fid, 'error')

    finally:
        remove_file_handlers(logger)
        remove_file_handlers(logger_warn)

    if table is not None and autoupload:
        api.upload_photometry(fid, delete_completed=True)
        api.set_photometry_status(fid, 'ingest')

    return table


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
    # group.add_argument('--no-imagematch', action='store_true', help="Disable ImageMatch.")
    group.add_argument('--rescale-static', action='store_true',
                       help="Rescale uncertainty using static FWHM=2.5x. Else use dynamic.")
    group.add_argument('--autoupload', action='store_true',
                       help="Automatically upload completed photometry to Flows website. "
                            "Only do this, if you know what you are doing!")
    group.add_argument('--fixposdiff', action='store_true',
                       help="Fix SN position during PSF photometry of difference image. "
                            "Useful when difference image is noisy.")
    group.add_argument('--rescale-off', action='store_false',
                        help='Turn off uncertainty rescaling.')
    group.add_argument('--wcstimeout', type=int, default=None, help="Timeout in Seconds for WCS.")
    args = parser.parse_args()

    # Ensure that all input has been given:
    if not args.fileid and not args.targetid and args.filter is None:
        parser.error("Please select either a specific FILEID .")

    # Setup logger:
    logging_level = parse_log_level(args)
    logger = create_logger(worker_name='parsing', log_level=logging_level)

    # Number of threads to use:
    threads = args.threads
    if threads <= 0:
        threads = multiprocessing.cpu_count()

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

    # Create function wrapper:
    process_fileid_wrapper = functools.partial(process_fileid, autoupload=args.autoupload,
                                               cm_timeout=args.wcstimeout, no_plots=args.noplots,
                                               rescale=args.rescale_off, rescale_dynamic=not args.rescale_static)

    if threads > 1:
        # process in parallel:
        with multiprocessing.Pool(threads) as pool:
            for result in tqdm.tqdm(pool.imap_unordered(process_fileid_wrapper, fileids), total=len(fileids)):
                if result is not None:
                    logger.info(f"finished:{result.meta['fileid']}")
                # We can do something with the partial results here (like calculate color terms!).
                pass

    else:
        # Only single thread so simply run it directly:
        for fid in fileids:
            # Setup logging:
            logger = create_logger(str(fid), log_level=logging_level)  # change to correct fid name for log.
            logger.info("=" * 72)
            logger.info(fid)
            logger.info("=" * 72)
            process_fileid_wrapper(fid)


if __name__ == '__main__':
    main()
