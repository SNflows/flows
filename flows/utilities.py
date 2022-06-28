#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Utility functions

.. codeauthor:: Rasmus Handberg <rasmush@phys.au.dk>
"""

import hashlib
import multiprocessing
import logging
import sys
from argparse import Namespace
from typing import Optional

from scipy.interpolate import interp2d


def get_filehash(fname):
    """Calculate SHA1-hash of file."""
    buf = 65536
    s = hashlib.sha1()
    with open(fname, 'rb') as fid:
        while True:
            data = fid.read(buf)
            if not data:
                break
            s.update(data)

    sha1sum = s.hexdigest().lower()
    if len(sha1sum) != 40:
        raise Exception("Invalid file hash")
    return sha1sum


def has_file_handler(logger):
    """Check if logger has one file handler."""
    return sum([type(l) is logging.FileHandler for l in logger.handlers]) > 0


def has_stream_handler(logger):
    """Check if logger has one stream handler."""
    return sum([type(l) is logging.StreamHandler for l in logger.handlers]) > 0


def remove_file_handlers(logger):
    """Remove file handler from logger."""
    for handler in logger.handlers:
        if type(handler) is logging.FileHandler:
            logger.removeHandler(handler)


def create_logger(worker_name=None, log_level: Optional[int] = None, log_file: str = None):
    formatter = logging.Formatter('[%(asctime)s| %(levelname)s| %(processName)s | %(module)s] %(message)s')
    if worker_name is not None:
        process = multiprocessing.current_process()
        process.name = worker_name

    logger = multiprocessing.get_logger()
    if log_level is not None:
        logger.setLevel(log_level)

    stream_handler = logging.StreamHandler(sys.stdout)
    stream_handler.setFormatter(formatter)

    # remove duplicated messages in the output
    if not has_stream_handler(logger):
        logger.addHandler(stream_handler)

    # check for log file and remove duplicated messages in the output
    if log_file is not None and not has_file_handler(logger):
        file_formatter = logging.Formatter(
            '%(asctime)s - %(levelname)s -%(module)s - %(message)s', "%Y-%m-%d %H:%M:%S")
        file_handler = logging.FileHandler(log_file, mode='w')
        file_handler.setFormatter(file_formatter)
        file_handler.setLevel(logging.INFO)
        logger.addHandler(file_handler)

    return logger


def create_warning_logger(log_file: str):
    """Create a logger for warnings."""
    logging.captureWarnings(True)
    logger_warn = logging.getLogger('py.warnings')

    file_formatter = logging.Formatter(
        '%(asctime)s - %(levelname)s -%(module)s - %(message)s', "%Y-%m-%d %H:%M:%S")
    file_handler = logging.FileHandler(log_file, mode='w')
    file_handler.setFormatter(file_formatter)
    file_handler.setLevel(logging.INFO)

    if not has_file_handler(logger_warn):
        logger_warn.addHandler(file_handler)

    return logger_warn


def parse_log_level(args: Namespace):
    logging_level = logging.INFO
    if args.quiet:
        logging_level = logging.WARNING
    elif args.debug:
        logging_level = logging.DEBUG
    return logging_level
