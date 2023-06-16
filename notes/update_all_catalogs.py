#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import logging
import sys
import os.path
import tqdm

if os.path.abspath('..') not in sys.path:
    sys.path.insert(0, os.path.abspath('..'))
import flows
from tendrils import api


# --------------------------------------------------------------------------------------------------
class TqdmLoggingHandler(logging.Handler):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def emit(self, record):
        try:
            msg = self.format(record)
            tqdm.tqdm.write(msg)
            self.flush()
        except (KeyboardInterrupt, SystemExit):  # pragma: no cover
            raise
        except:  # noqa: E722, pragma: no cover
            self.handleError(record)


# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    # Setup logging:
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    console = TqdmLoggingHandler()
    console.setFormatter(formatter)
    logger = logging.getLogger('flows')
    if not logger.hasHandlers():
        logger.addHandler(console)
    logger.setLevel(logging.INFO)

    # Do it by status, just to prioritize things a bit:
    for tgtstatus in ('target', 'candidate', 'rejected'):
        targetids = sorted([tgt['targetid'] for tgt in api.get_targets() if tgt['target_status'] == tgtstatus])[
                    ::-1]

        for targetid in tqdm.tqdm(targetids, desc=tgtstatus):
            donefile = f"catalog_updates/{targetid:05d}.done"
            if not os.path.exists(donefile):
                try:
                    flows.catalogs.download_catalog(targetid, update_existing=True)
                except:
                    logger.exception("targetid=%d", targetid)
                else:
                    open(donefile, 'w').close()
