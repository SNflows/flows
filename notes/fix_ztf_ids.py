#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import os.path
from tqdm import tqdm
from astropy.coordinates import SkyCoord

if os.path.abspath('..') not in sys.path:
    sys.path.insert(0, os.path.abspath('..'))
import flows
from tendrils import api

# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    with flows.aadc_db.AADC_DB() as db:
        for target in tqdm(api.get_targets()):
            if target['target_status'] == 'rejected':
                continue

            targetid = target['targetid']
            coord = SkyCoord(ra=target['ra'], dec=target['decl'], unit='deg', frame='icrs')
            dd = target['discovery_date']

            # Query for the ZTF id:
            ztf_id = flows.ztf.query_ztf_id(coord, discovery_date=dd)

            # If the ZTF id is not the same as we have currently, update it in the database:
            if ztf_id != target['ztf_id']:
                print(target)
                print(ztf_id)
                print("******* NEEDS UPDATE ******")

                db.cursor.execute("UPDATE flows.targets SET ztf_id=%s WHERE targetid=%s;", (ztf_id, targetid))
                db.conn.commit()
