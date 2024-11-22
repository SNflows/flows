#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import logging
import os
import datetime
import shutil
import glob
import re
import requests
import urllib3
import ssl
import sys
from tqdm import tqdm
if os.path.abspath(os.path.join(os.path.dirname(__file__), 'flows')) not in sys.path:
    sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), 'flows')))
import flows
from flows.target import Target
from flows.aadc_db import AADC_DB
from tendrils import api

#---------
# https://stackoverflow.com/questions/66689696/urllib3-error-ssl-wrong-signature-type
class SslOldHttpAdapter(requests.adapters.HTTPAdapter):
    def init_poolmanager(self, connections, maxsize, block=False):
        ctx = ssl.create_default_context()
        ctx.set_ciphers('DEFAULT@SECLEVEL=1')
        self.poolmanager = urllib3.poolmanager.PoolManager(
            ssl_version=ssl.PROTOCOL_TLS,
            ssl_context=ctx)

#---------
if __name__ == '__main__':

    # Parse command line arguments:
    parser = argparse.ArgumentParser(description='Download ZTF photometry.')
    parser.add_argument('-d', '--debug', help='Print debug messages.', action='store_true')
    parser.add_argument('-q', '--quiet', help='Only report warnings and errors.', action='store_true')
    parser.add_argument('--fast', action='store_true')
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

    today = datetime.datetime.utcnow()
    tomorrow = datetime.datetime.utcnow() + datetime.timedelta(hours=24)
    yesterday = datetime.datetime.utcnow() - datetime.timedelta(hours=24)

    # Automatically reject targets that are more than 2 months old:



    with AADC_DB() as db:
        db.cursor.execute("""UPDATE flows.targets SET target_status='rejected' WHERE
			target_status='candidate'
			AND inserted < (NOW() at time zone 'utc') - '3 months'::interval;""") # TODO: Change to last_modified!
        n_rejected = db.cursor.rowcount
        db.conn.commit()
        logger.info("%d targets automatically rejected.", n_rejected)

    # Sort the targets by when they were inserted:
    targets = api.get_targets()
    targets = sorted(targets, key=lambda x: x['inserted'])[::-1]

    tqdm_settings = {'disable': None if logger.isEnabledFor(logging.INFO) else True}

    # Start requests session and allow older SSL settings
    # since sdss.org is not that modern:
    ses = requests.Session()
    ses.mount('https://skyserver.sdss.org', SslOldHttpAdapter())

    regex_year = re.compile(r'^(\d+)')
    for tgt in tqdm(targets, **tqdm_settings):

        m = regex_year.match(tgt['target_name'])
        year = int(m.group(1))

        outdir = os.path.join('/archive/candidates', str(year), tgt['target_name'])

        # If running in fast mode, jump over everything
        if args.fast and (tgt['inserted'] < yesterday or os.path.isdir(outdir)):
            continue

        if tgt['target_status'] not in ('candidate', 'target'):
            # Cleanup any leftover directories:
            if os.path.isdir(outdir):
                logger.info("Deleting old directory for %s...", tgt['target_name'])
                shutil.rmtree(outdir)
            continue

        logger.info("Creating plots for %s...", tgt['target_name'])
        os.makedirs(outdir, exist_ok=True)

        tgtobj = Target.from_tid(tgt['targetid'])

        # FIXME: Limiting to only siteid=2 (VLT) right now!
        flows.visibility(target=tgtobj, siteid=2, date=today.strftime('%Y%m%d'), output=outdir, overwrite=False)

        # Create plots for tomorrow already:
        flows.visibility(target=tgtobj, siteid=2, date=tomorrow.strftime('%Y%m%d'), output=outdir, overwrite=False)

        regex_fname = re.compile(r'^visibility_%s_(\d+)_site\d+\.png$' % tgt['target_name'])
        for f in glob.iglob(os.path.join(outdir, 'visibility_*.png')):
            bname = os.path.basename(f)

            m = regex_fname.match(bname)
            if m:
                file_time = datetime.datetime.strptime(m.group(1), '%Y%m%d')
                file_age = file_time - today

                if file_age < datetime.timedelta(hours=-48):
                    logger.info("Deleting old file: %s", bname)
                    os.remove(f)

            else:
                logger.warning("Deleting non-matching file: %s", f)
                os.remove(f)

        #------------------------------------------------------------------------------------------

        filename = tgt['target_name'] + '_sdss_thumb.jpg'
        filepath = os.path.join(outdir, filename)
        if not os.path.isfile(filepath):
            logger.info("Downloading SDSS thumbnail for %s...", tgt['target_name'])
            params = {
                'ra': tgt['ra'],
                'dec': tgt['decl'],
                'scale': 0.8, # 1.79224,
                'width': 150,
                'height': 150,
                'opt': 'GI'
            }

            r = ses.get('https://skyserver.sdss.org/dr16/SkyServerWS/ImgCutout/getjpeg',
                        params=params, stream=True)
            if r.status_code == 404: # This is allowed!
                if os.path.isfile(filepath):
                    os.remove(filepath)
            else:
                r.raise_for_status()

                # Save the image to disk:
                with open(filepath, 'wb') as fd:
                    for chunk in r.iter_content(chunk_size=128):
                        fd.write(chunk)
