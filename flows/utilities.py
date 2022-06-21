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


def create_logger(worker_name=None):
    if worker_name is not None:
        process = multiprocessing.current_process()
        process.name = worker_name

    logger = multiprocessing.get_logger()
    logger.setLevel(logging.INFO)
    formatter = logging.Formatter('[%(asctime)s| %(levelname)s| %(processName)s | %(module)s] %(message)s')
    handler = logging.StreamHandler(sys.stdout)
    handler.setFormatter(formatter)

    # remove duplicated messages in the output
    if not len(logger.handlers):
        logger.addHandler(handler)


    #    logger.level = level_name


    return logger
