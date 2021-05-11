#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from tqdm import tqdm
from astropy.coordinates import SkyCoord
import flows

#--------------------------------------------------------------------------------------------------
if __name__ == '__main__':
	with flows.aadc_db.AADC_DB() as db:
		for target in tqdm(flows.api.get_targets()):
			if target['target_status'] == 'rejected':
				continue

			targetid = target['targetid']
			coord = SkyCoord(ra=target['ra'], dec=target['decl'], unit='deg', frame='icrs')
			ztf_id = flows.ztf.query_ztf_id(coord)

			if ztf_id != target['ztf_id']:
				print(target)
				print(ztf_id)
				print("******* NEEDS UPDATE ******")

				db.cursor.execute("UPDATE flows.targets SET ztf_id=%s WHERE targetid=%s;", (ztf_id, targetid))
				db.conn.commit()
