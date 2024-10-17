#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os.path
import getpass
import numpy as np
from astropy.table import Table
from astropy.units import Quantity
import sys
if sys.path[0] != os.path.abspath(os.path.join(os.path.dirname(__file__), 'flows')):
    sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), 'flows')))
from flows.aadc_db import AADC_DB
from tendrils.utils import load_config

if __name__ == '__main__':
    with AADC_DB() as db:

        #db.cursor.execute("SELECT fileid,path,targetid FROM flows.files LEFT JOIN flows.photometry_details d ON d.fileid_phot=files.fileid WHERE files.datatype=2 AND d.fileid_phot IS NULL;")
        #
        #for row in db.cursor.fetchall():
        #
        #	fileid_phot = row['fileid']
        #	filepath = os.path.join('/archive/photometry/', row['path'])
        #
        #	tab = Table.read(filepath, format='ascii.ecsv')
        #
        #	indx_raw = (tab['starid'] == 0)
        #	indx_sub = (tab['starid'] == -1)
        #	indx_ref = (tab['starid'] > 0)
        #
        #	phot_details = {
        #		'fileid_phot': fileid_phot,
        #		'fileid_img': int(tab.meta['fileid']),
        #		'fileid_template': tab.meta['template'],
        #		'fileid_diffimg': None if 'diffimg' not in tab.meta else tab.meta['diffimg'],
        #		'obstime': tab.meta['obstime-bmjd'],
        #		'mag_raw': float(tab[indx_raw]['mag']),
        #		'mag_raw_error': float(tab[indx_raw]['mag_error']),
        #		'mag_sub': None if not any(indx_sub) else float(tab[indx_sub]['mag']),
        #		'mag_sub_error': None if not any(indx_sub) else float(tab[indx_sub]['mag_error']),
        #		'zeropoint': float(tab.meta['zp']),
        #		'zeropoint_error': None if 'zp_error' not in tab.meta else float(tab.meta['zp_error']),
        #		'zeropoint_diff': None if 'zp_diff' not in tab.meta else float(tab.meta['zp_diff']),
        #		'fwhm': None if 'fwhm' not in tab.meta else tab.meta['fwhm'],
        #		'seeing': None if 'seeing' not in tab.meta else tab.meta['seeing'],
        #		'references_detected': int(np.sum(indx_ref)),
        #		'used_for_epsf': int(np.sum(tab['used_for_epsf'])),
        #		'faintest_reference_detected': float(np.max(tab[indx_ref]['mag'])),
        #		'pipeline_version': tab.meta['version'],
        #	}
        #
        #	for key in ('fwhm', 'seeing'):
        #		if isinstance(phot_details[key], Quantity):
        #			phot_details[key] = float(phot_details[key].value)
        #
        #	print(phot_details)
        #	print(row['targetid'])
        #
        #	db.cursor.execute("""INSERT INTO flows.photometry_details (
        #		fileid_phot,
        #		fileid_img,
        #		fileid_template,
        #		fileid_diffimg,
        #		obstime_bmjd,
        #		mag_raw,
        #		mag_raw_error,
        #		mag_sub,
        #		mag_sub_error,
        #		zeropoint,
        #		zeropoint_error,
        #		zeropoint_diff,
        #		fwhm,
        #		seeing,
        #		references_detected,
        #		used_for_epsf,
        #		faintest_reference_detected,
        #		pipeline_version
        #	) VALUES (
        #		%(fileid_phot)s,
        #		%(fileid_img)s,
        #		%(fileid_template)s,
        #		%(fileid_diffimg)s,
        #		%(obstime)s,
        #		%(mag_raw)s,
        #		%(mag_raw_error)s,
        #		%(mag_sub)s,
        #		%(mag_sub_error)s,
        #		%(zeropoint)s,
        #		%(zeropoint_error)s,
        #		%(zeropoint_diff)s,
        #		%(fwhm)s,
        #		%(seeing)s,
        #		%(references_detected)s,
        #		%(used_for_epsf)s,
        #		%(faintest_reference_detected)s,
        #		%(pipeline_version)s
        #	);""", phot_details)
        #
        #db.conn.commit()

        db.cursor.execute("SELECT fileid,path,targetid FROM flows.files LEFT JOIN flows.photometry_details d ON d.fileid_phot=files.fileid WHERE files.datatype=2 AND d.faintest_reference_detected='NaN'::real;")
        for row in db.cursor.fetchall():

            fileid_phot = row['fileid']
            filepath = os.path.join('/archive/photometry/', row['path'])

            tab = Table.read(filepath, format='ascii.ecsv')
            print(len(tab))

            indx_ref = (tab['starid'] > 0)

            frd = float(np.nanmax(tab[indx_ref]['mag']))
            if np.isnan(frd):
                frd = None

            print(fileid_phot, frd)

            db.cursor.execute("UPDATE flows.photometry_details SET faintest_reference_detected=%s WHERE fileid_phot=%s;", [frd, fileid_phot])

        db.conn.commit()
