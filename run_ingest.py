#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Run ingest of data files in Flows system.

This code is obviously only meant to run on the central
Flows systems, and will not work outside of that environment.
@TODO: Refactor Database out of this script.
.. codeauthor:: Rasmus Handberg <rasmush@phys.au.dk>
"""

import logging
import sys
import os.path
import glob
import subprocess
import warnings
import numpy as np
import astropy.coordinates as coords
from astropy.table import Table
import shutil
import gzip
from zipfile import ZipFile
import tempfile
import re
from collections import defaultdict
if os.path.abspath(os.path.join(os.path.dirname(__file__), 'flows')) not in sys.path:
    sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), 'flows')))
from flows.aadc_db import AADC_DB
from flows.plots import plt, plot_image
from flows.load_image import load_image
from flows.utilities import get_filehash


# --------------------------------------------------------------------------------------------------
def flows_get_archive_from_path(fname, archives_list=None):
    """
    Translate full path into AADC archive identifier and relative path.

    It is highly recommended to provide the list with that call
    to this function since it will involve a query to the database
    at every call.
    """

    archive = None
    relpath = None

    # Get list of archives, if not provided with call:
    if archives_list is None:
        with AADC_DB() as db:
            db.cursor.execute("SELECT archive,path FROM aadc.files_archives;")
            archives_list = db.cursor.fetchall()

    # Make sure folder is absolute path
    folder = os.path.abspath(fname)

    # Loop through the defined archives and find one that matches:
    for opt in archives_list:
        archive_path = opt['path']
        if archive_path is not None and archive_path != '':
            archive_path = archive_path.rstrip('/\\') + os.path.sep
            if folder.startswith(archive_path):
                archive = int(opt['archive'])
                relpath = folder[len(archive_path):].replace('\\', '/')
                break

    # We did not find anything:
    if archive is None:
        raise RuntimeError("File not in registred archive")

    return archive, relpath


# --------------------------------------------------------------------------------------------------
def optipng(fpath):
    subprocess.run(['optipng', '-preserve', '-quiet', fpath])

# --------------------------------------------------------------------------------------------------
class CounterFilter(logging.Filter):
    """
    A logging filter which counts the number of log records in each level.
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.counter = defaultdict(int)

    def filter(self, record):  # noqa: A003
        self.counter[record.levelname] += 1
        return True


# --------------------------------------------------------------------------------------------------
def create_plot(filepath, target_coord=None, target_position=None):
    output_fpath = os.path.abspath(re.sub(r'\.fits(\.gz)?$', '', filepath) + '.png')

    img = load_image(filepath, target_coord=target_coord)

    fig = plt.figure(figsize=(12, 12))
    ax = fig.add_subplot(111)
    plot_image(img.clean, ax=ax, scale='linear', percentile=[5, 99], cbar='right')
    if target_position is not None:
        ax.scatter(target_position[0], target_position[1], marker='+', s=20, c='r', label='Target')
    fig.savefig(output_fpath, bbox_inches='tight')
    plt.close(fig)

    optipng(output_fpath)


# --------------------------------------------------------------------------------------------------
def ingest_from_inbox():
    rootdir_inbox = '/data/flows/archive/inbox'
    rootdir = '/data/flows/archive/archive'

    logger = logging.getLogger(__name__)

    # Check that root directories are available:
    if not os.path.isdir(rootdir_inbox):
        raise FileNotFoundError("INBOX does not exists")
    if not os.path.isdir(rootdir):
        raise FileNotFoundError("ARCHIVE does not exists")

    with AADC_DB() as db:
        # Get list of archives:
        db.cursor.execute("SELECT archive,path FROM aadc.files_archives;")
        archives_list = db.cursor.fetchall()

        # Get list of all available filters:
        db.cursor.execute("SELECT photfilter FROM flows.photfilters;")
        all_filters = set([row['photfilter'] for row in db.cursor.fetchall()])

        for inputtype in ('science', 'templates', 'subtracted', 'replace'):  #
            for fpath in glob.iglob(os.path.join(rootdir_inbox, '*', inputtype, '*')):
                logger.info("=" * 72)
                logger.info(fpath)

                # Find the uploadlog corresponding to this file:
                db.cursor.execute("SELECT logid,dont_ingest_again FROM flows.uploadlog WHERE uploadpath=%s;", [os.path.relpath(fpath, rootdir_inbox)])
                row = db.cursor.fetchone()
                if row is not None:
                    uploadlogid = row['logid']
                    if row['dont_ingest_again']:
                        logger.info("Skipping file")
                        continue
                else:
                    uploadlogid = None
                logger.info("Uploadlog ID: %s", uploadlogid)

                # Only accept FITS file, or already compressed FITS files:
                if not fpath.endswith('.fits') and not fpath.endswith('.fits.gz'):
                    if uploadlogid:
                        db.cursor.execute("UPDATE flows.uploadlog SET status='Invalid file type',dont_ingest_again=TRUE WHERE logid=%s;", [uploadlogid])
                        db.conn.commit()
                    logger.error("Invalid file type: %s", os.path.relpath(fpath, rootdir_inbox))
                    continue

                # Get the name of the directory:
                # Not pretty, but it works...
                target_dirname = fpath[len(rootdir_inbox) + 1:]
                target_dirname = target_dirname.split(os.path.sep)[0]

                # Convert directory name to target
                db.cursor.execute("SELECT targetid,target_name,ra,decl FROM flows.targets WHERE target_name=%s;",
                                  [target_dirname])
                row = db.cursor.fetchone()
                if row is None:
                    logger.error('Could not find target: %s', target_dirname)
                    continue
                targetid = row['targetid']
                targetname = row['target_name']
                target_radec = [[row['ra'], row['decl']]]
                target_coord = coords.SkyCoord(ra=row['ra'], dec=row['decl'], unit='deg', frame='icrs')

                if not fpath.endswith('.gz'):
                    # Gzip the FITS file:
                    with open(fpath, 'rb') as f_in:
                        with gzip.open(fpath + '.gz', 'wb') as f_out:
                            f_out.writelines(f_in)

                    # We should now have a Gzip file instead:
                    if os.path.isfile(fpath) and os.path.isfile(fpath + '.gz'):
                        # Update the log of this file:
                        if uploadlogid:
                            db.cursor.execute("UPDATE flows.uploadlog SET uploadpath=%s WHERE logid=%s;",
                                              [os.path.relpath(fpath + '.gz', rootdir_inbox), uploadlogid])
                            db.conn.commit()

                        os.remove(fpath)
                        fpath += '.gz'
                    else:
                        raise RuntimeError("Gzip file was not created correctly")

                version = 1
                if inputtype == 'science':
                    newpath = os.path.join(rootdir, targetname, os.path.basename(fpath))
                    datatype = 1
                elif inputtype == 'templates':
                    newpath = os.path.join(rootdir, targetname, inputtype, os.path.basename(fpath))
                    datatype = 3
                elif inputtype == 'subtracted':
                    newpath = os.path.join(rootdir, targetname, inputtype, os.path.basename(fpath))
                    datatype = 4

                    original_fname = os.path.basename(fpath).replace('diff.fits', '.fits')
                    db.cursor.execute(
                        "SELECT fileid FROM flows.files WHERE targetid=%s AND datatype=1 AND path LIKE %s;",
                        [targetid, '%/' + original_fname])
                    subtracted_original_fileid = db.cursor.fetchone()
                    if subtracted_original_fileid is None:
                        if uploadlogid:
                            db.cursor.execute(
                                "UPDATE flows.uploadlog SET status='original science image not found' WHERE logid=%s;",
                                [uploadlogid])
                            db.conn.commit()
                        logger.error("ORIGINAL SCIENCE IMAGE COULD NOT BE FOUND: %s", os.path.basename(fpath))
                        continue
                    else:
                        subtracted_original_fileid = subtracted_original_fileid[0]

                elif inputtype == 'replace':
                    bname = os.path.basename(fpath)
                    m = re.match(r'^(\d+)_v(\d+)\.fits(\.gz)?$', bname)
                    if m:
                        replaceid = int(m.group(1))
                        version = int(m.group(2))

                        db.cursor.execute("SELECT datatype,path,version FROM flows.files WHERE fileid=%s;", [replaceid])
                        row = db.cursor.fetchone()
                        if row is None:
                            logger.error("Unknown fileid to be replaced: %s", bname)
                            continue
                        datatype = row['datatype']
                        subdir = {1: '', 4: 'subtracted'}[datatype]

                        if version != row['version'] + 1:
                            logger.error("Mismatch in versions: old=%d, new=%d", row['version'], version)
                            continue

                        newfilename = re.sub(r'(_v\d+)?\.fits(\.gz)?$', r'_v{version:d}.fits\2'.format(version=version),
                                             os.path.basename(row['path']))
                        newpath = os.path.join(rootdir, targetname, subdir, newfilename)

                        if datatype == 4:
                            db.cursor.execute(
                                "SELECT associd FROM flows.files_cross_assoc INNER JOIN flows.files ON files.fileid=files_cross_assoc.associd WHERE files_cross_assoc.fileid=%s AND datatype=1;",
                                [replaceid])
                            subtracted_original_fileid = db.cursor.fetchone()
                            if subtracted_original_fileid is None:
                                if uploadlogid:
                                    db.cursor.execute(
                                        "UPDATE flows.uploadlog SET status='original science image not found' WHERE logid=%s;",
                                        [uploadlogid])
                                    db.conn.commit()
                                logger.error("ORIGINAL SCIENCE IMAGE COULD NOT BE FOUND: %s", os.path.basename(fpath))
                                continue
                            else:
                                subtracted_original_fileid = subtracted_original_fileid[0]
                    else:
                        logger.error("Invalid replace file name: %s", bname)
                        continue

                else:
                    raise RuntimeError("Not understood, Captain")

                logger.info(newpath)

                if os.path.exists(newpath):
                    logger.error("Already exists")
                    if uploadlogid:
                        db.cursor.execute(
                            "UPDATE flows.uploadlog SET status='Already exists: file name' WHERE logid=%s;",
                            [uploadlogid])
                        db.conn.commit()
                    continue

                archive, relpath = flows_get_archive_from_path(newpath, archives_list)

                db.cursor.execute("SELECT fileid FROM flows.files WHERE archive=%s AND path=%s;", [archive, relpath])
                if db.cursor.fetchone() is not None:
                    logger.error("ALREADY DONE")
                    continue

                # Calculate filehash of the file being stored:
                filehash = get_filehash(fpath)

                # Check that the file does not already exist:
                db.cursor.execute("SELECT fileid FROM flows.files WHERE filehash=%s;", [filehash])
                if db.cursor.fetchone() is not None:
                    logger.error("ALREADY DONE: Filehash")
                    if uploadlogid:
                        db.cursor.execute(
                            "UPDATE flows.uploadlog SET status='Already exists: filehash' WHERE logid=%s;",
                            [uploadlogid])
                        db.conn.commit()
                    continue

                # Try to load the image using the same function as the pipeline would:
                try:
                    img = load_image(fpath, target_coord=target_coord)
                except Exception as e:  # pragma: no cover
                    logger.exception("Could not load FITS image")
                    if uploadlogid:
                        errmsg = str(e) if not hasattr(e, 'message') else str(e.message)
                        db.cursor.execute("UPDATE flows.uploadlog SET status=%s WHERE logid=%s;",
                                          ['Load Image Error: ' + errmsg, uploadlogid])
                        db.conn.commit()
                    continue

                # Use the WCS in the file to calculate the pixel-positon of the target:
                with warnings.catch_warnings():
                    warnings.filterwarnings('ignore', category=RuntimeWarning, message='All-NaN slice encountered')
                try:
                    target_pixels = img.wcs.all_world2pix(target_radec, 0).flatten()
                except:  # noqa: E722, pragma: no cover
                    logger.exception("Could not find target position using the WCS.")
                    if uploadlogid:
                        errmsg = "Could not find target position using the WCS."
                        db.cursor.execute("UPDATE flows.uploadlog SET status=%s,dont_ingest_again=TRUE WHERE logid=%s;", [errmsg, uploadlogid])
                        db.conn.commit()
                    continue

                # Check that the position of the target actually falls within
                # the pixels of the image:
                logger.debug("image.shape = %s", img.shape)
                if len(img.shape) != 2:
                    errmsg = "Invalid image dimensions."
                    logger.error(errmsg)
                    if uploadlogid:
                        db.cursor.execute("UPDATE flows.uploadlog SET status=%s,dont_ingest_again=TRUE WHERE logid=%s;", [errmsg, uploadlogid])
                        db.conn.commit()
                    continue

                logger.debug("target_pixels = %s", target_pixels)
                if len(target_pixels) != 2 or np.any(~np.isfinite(target_pixels)):
                    errmsg = "Invalid extracted target position. Check the WCS."
                    logger.error(errmsg)
                    if uploadlogid:
                        db.cursor.execute("UPDATE flows.uploadlog SET status=%s,dont_ingest_again=TRUE WHERE logid=%s;", [errmsg, uploadlogid])
                        db.conn.commit()
                    continue

                if target_pixels[0] < -0.5 or target_pixels[1] < -0.5 \
                    or target_pixels[0] > img.shape[1]-0.5 or target_pixels[1] > img.shape[0]-0.5:
                    logger.error("Target position does not fall within image. Check the WCS.")
                    if uploadlogid:
                        errmsg = "Target position does not fall within image. Check the WCS."
                        db.cursor.execute("UPDATE flows.uploadlog SET status=%s,dont_ingest_again=TRUE WHERE logid=%s;", [errmsg, uploadlogid])
                        db.conn.commit()
                    continue

                # Check that the site was found:
                if img.site is None or img.site['siteid'] is None:
                    logger.error("Unknown SITE")
                    if uploadlogid:
                        errmsg = "Unknown site"
                        db.cursor.execute("UPDATE flows.uploadlog SET status=%s WHERE logid=%s;", [errmsg, uploadlogid])
                        db.conn.commit()
                    continue

                # Check that the extracted photometric filter is valid:
                if img.photfilter not in all_filters:
                    logger.error("Unknown PHOTFILTER: %s", img.photfilter)
                    if uploadlogid:
                        errmsg = "Unknown PHOTFILTER: '" + str(img.photfilter) + "'"
                        db.cursor.execute("UPDATE flows.uploadlog SET status=%s WHERE logid=%s;", [errmsg, uploadlogid])
                        db.conn.commit()
                    continue

                # Do a deep check to ensure that there is not already another file with the same
                # properties (target, datatype, site, filter) taken at the same time:
                # TODO: Look at the actual overlap with the database, instead of just overlap
                #       with the central value. This way of doing it is more forgiving.
                obstime = img.obstime.utc.mjd
                if inputtype != 'replace':
                    db.cursor.execute("SELECT fileid FROM flows.files WHERE targetid=%s AND datatype=%s AND site=%s AND photfilter=%s AND obstime BETWEEN %s AND %s;", [
                        targetid,
                        datatype,
                        img.site['siteid'],
                        img.photfilter,
                        float(obstime - 0.5 * img.exptime/86400),
                        float(obstime + 0.5 * img.exptime / 86400), ])
                    if db.cursor.fetchone() is not None:
                        logger.error("ALREADY DONE: Deep check")
                        if uploadlogid:
                            db.cursor.execute(
                                "UPDATE flows.uploadlog SET status='Already exists: deep check' WHERE logid=%s;",
                                [uploadlogid])
                            db.conn.commit()
                        continue

                try:
                    # Copy the file to its new home:
                    os.makedirs(os.path.dirname(newpath), exist_ok=True)
                    shutil.copy(fpath, newpath)

                    # Set file and directory permissions:
                    # TODO: Can this not be handled in a more elegant way?
                    os.chmod(os.path.dirname(newpath), 0o2750)
                    os.chmod(newpath, 0o0440)

                    filesize = os.path.getsize(fpath)

                    if not fpath.endswith('-e00.fits'):
                        create_plot(newpath, target_coord=target_coord, target_position=target_pixels)

                    db.cursor.execute("INSERT INTO flows.files (archive,path,targetid,datatype,site,filesize,filehash,obstime,photfilter,exptime,version,available) VALUES (%(archive)s,%(relpath)s,%(targetid)s,%(datatype)s,%(site)s,%(filesize)s,%(filehash)s,%(obstime)s,%(photfilter)s,%(exptime)s,%(version)s,1) RETURNING fileid;", {
                        'archive': archive,
                        'relpath': relpath,
                        'targetid': targetid,
                        'datatype': datatype,
                         'site': img.site['siteid'], 'filesize': filesize, 'filehash': filehash, 'obstime': obstime,
                         'photfilter': img.photfilter, 'exptime': img.exptime, 'version': version})
                    fileid = db.cursor.fetchone()[0]

                    if datatype == 4:
                        db.cursor.execute("INSERT INTO flows.files_cross_assoc (fileid,associd) VALUES (%s,%s);",
                                          [fileid, subtracted_original_fileid])

                    if inputtype == 'replace':
                        db.cursor.execute("UPDATE flows.files SET newest_version=FALSE WHERE fileid=%s;", [replaceid])

                    if uploadlogid:
                        db.cursor.execute("UPDATE flows.uploadlog SET fileid=%s,status='ok' WHERE logid=%s;",
                                          [fileid, uploadlogid])

                    db.conn.commit()

                except (KeyboardInterrupt, SystemExit):
                    db.conn.rollback()
                    if os.path.exists(newpath):
                        os.remove(newpath)
                    logger.warning("Stopped by user or system")
                    return

                except Exception as e: # noqa: E722, pragma: no cover
                    db.conn.rollback()
                    if os.path.exists(newpath):
                        os.remove(newpath)
                    logger.exception("%s: Could not ingest file", fpath)
                    if uploadlogid:
                        msg = str(e) if not hasattr(e, 'message') else str(e.message)
                        db.cursor.execute("UPDATE flows.uploadlog SET status=%s WHERE logid=%s;", [msg, uploadlogid])
                        db.conn.commit()

                else:
                    logger.info("DELETE THE ORIGINAL FILE")
                    if os.path.isfile(newpath):
                        os.remove(fpath)
                    if uploadlogid:
                        db.cursor.execute("UPDATE flows.uploadlog SET uploadpath=NULL WHERE logid=%s;", [uploadlogid])
                        db.conn.commit()


# --------------------------------------------------------------------------------------------------
def ingest_photometry_from_inbox():
    rootdir_inbox = '/data/flows/archive/inbox'
    rootdir_archive = '/data/flows/archive/archive_photometry'

    logger = logging.getLogger(__name__)

    # Check that root directories are available:
    if not os.path.isdir(rootdir_inbox):
        raise FileNotFoundError("INBOX does not exists")
    if not os.path.isdir(rootdir_archive):
        raise FileNotFoundError("ARCHIVE does not exists")

    with AADC_DB() as db:
        # Get list of archives:
        db.cursor.execute("SELECT archive,path FROM aadc.files_archives;")
        archives_list = db.cursor.fetchall()

        for fpath in glob.iglob(os.path.join(rootdir_inbox, '*', 'photometry', '*')):
            logger.info("=" * 72)
            logger.info(fpath)

            # Find the uploadlog corresponding to this file:
            db.cursor.execute("SELECT logid FROM flows.uploadlog WHERE uploadpath=%s;",
                              [os.path.relpath(fpath, rootdir_inbox)])
            row = db.cursor.fetchone()
            if row is not None:
                uploadlogid = row['logid']
            else:
                uploadlogid = None

            # Only accept FITS file, or already compressed FITS files:
            if not fpath.endswith('.zip'):
                if uploadlogid:
                    db.cursor.execute("UPDATE flows.uploadlog SET status='Invalid file type' WHERE logid=%s;",
                                      [uploadlogid])
                    db.conn.commit()
                logger.error("Invalid file type: %s", os.path.relpath(fpath, rootdir_inbox))
                continue

            # Get the name of the directory:
            # Not pretty, but it works...
            target_dirname = fpath[len(rootdir_inbox) + 1:]
            target_dirname = target_dirname.split(os.path.sep)[0]

            # Convert directory name to target
            db.cursor.execute("SELECT targetid,target_name FROM flows.targets WHERE target_name=%s;", [target_dirname])
            row = db.cursor.fetchone()
            if row is None:
                logger.error('Could not find target: %s', target_dirname)
                continue
            targetid = row['targetid']
            targetname = row['target_name']

            newpath = None
            try:
                with tempfile.TemporaryDirectory() as tmpdir:
                    #
                    tmpphotfile = os.path.join(tmpdir, 'photometry.ecsv')

                    # Extract the ZIP file:
                    with ZipFile(fpath, mode='r') as myzip:
                        for member in myzip.infolist():
                            # Remove any directory structure from the zip file:
                            if member.filename.endswith('/'):  # member.is_dir()
                                continue
                            member.filename = os.path.basename(member.filename)

                            # Due to security considerations, we only allow specific files
                            # to be extracted:
                            if member.filename == 'photometry.ecsv':
                                myzip.extract(member, path=tmpdir)
                            elif member.filename.endswith('.png') or member.filename.endswith('.log'):
                                myzip.extract(member, path=tmpdir)

                    # Check that the photometry ECSV file at least exists:
                    if not os.path.isfile(tmpphotfile):
                        raise FileNotFoundError("Photometry is not found")

                    # Load photometry table:
                    tab = Table.read(tmpphotfile, format='ascii.ecsv')
                    fileid_img = int(tab.meta['fileid'])
                    targetid_table = int(tab.meta['targetid'])

                    assert targetid_table == targetid

                    # Find out which version number to assign to file:
                    db.cursor.execute(
                        "SELECT MAX(files.version) AS latest_version FROM flows.files_cross_assoc fca INNER JOIN flows.files ON fca.fileid=files.fileid WHERE fca.associd=%s AND files.datatype=2;",
                        [fileid_img, ])
                    latest_version = db.cursor.fetchone()
                    if latest_version[0] is None:
                        new_version = 1
                    else:
                        new_version = latest_version[0] + 1

                    # Create a new path and filename that is slightly more descriptive:
                    newpath = os.path.join(rootdir_archive, targetname, f'{fileid_img:05d}', f'v{new_version:02d}',
                                           f'photometry-{targetname:s}-{fileid_img:05d}-v{new_version:02d}.ecsv')
                    logger.info(newpath)

                    if os.path.exists(newpath):
                        logger.error("Already exists")
                        if uploadlogid:
                            db.cursor.execute(
                                "UPDATE flows.uploadlog SET status='Already exists: file name' WHERE logid=%s;",
                                [uploadlogid])
                            db.conn.commit()
                        continue

                    archive, relpath = flows_get_archive_from_path(newpath, archives_list)

                    db.cursor.execute("SELECT fileid FROM flows.files WHERE archive=%s AND path=%s;",
                                      [archive, relpath])
                    if db.cursor.fetchone() is not None:
                        logger.error("ALREADY DONE")
                        continue

                    db.cursor.execute("SELECT * FROM flows.files WHERE fileid=%s;", [fileid_img])
                    row = db.cursor.fetchone()
                    site = row['site']

                    assert targetid == row['targetid']
                    assert tab.meta['photfilter'] == row['photfilter']

                    # Optimize all the PNG files in the temp directory:
                    for f in glob.iglob(os.path.join(tmpdir, '*.png')):
                        logger.debug("Running optipng on %s", f)
                        optipng(f)

                    # Copy the full directory to its new home:
                    shutil.copytree(tmpdir, os.path.dirname(newpath))
                    os.rename(os.path.join(os.path.dirname(newpath), 'photometry.ecsv'), newpath)

                # Get information about file:
                filesize = os.path.getsize(newpath)
                filehash = get_filehash(newpath)

                db.cursor.execute("INSERT INTO flows.files (archive,path,targetid,datatype,site,filesize,filehash,obstime,photfilter,version,available) VALUES (%(archive)s,%(relpath)s,%(targetid)s,%(datatype)s,%(site)s,%(filesize)s,%(filehash)s,%(obstime)s,%(photfilter)s,%(version)s,1) RETURNING fileid;", {
                    'archive': archive,
                    'relpath': relpath,
                    'targetid': targetid,
                    'datatype': 2,
                    'site': site,
                     'filesize': filesize, 'filehash': filehash, 'obstime': tab.meta['obstime-bmjd'],
                     'photfilter': tab.meta['photfilter'], 'version': new_version})
                fileid = db.cursor.fetchone()[0]

                # Add dependencies:
                db.cursor.execute("INSERT INTO flows.files_cross_assoc (fileid,associd) VALUES (%s,%s);",
                                  [fileid, fileid_img])
                if tab.meta['template'] is not None:
                    db.cursor.execute("INSERT INTO flows.files_cross_assoc (fileid,associd) VALUES (%s,%s);",
                                      [fileid, tab.meta['template']])
                if tab.meta['diffimg'] is not None:
                    db.cursor.execute("INSERT INTO flows.files_cross_assoc (fileid,associd) VALUES (%s,%s);",
                                      [fileid, tab.meta['diffimg']])

                indx_raw = (tab['starid'] == 0)
                indx_sub = (tab['starid'] == -1)
                indx_ref = (tab['starid'] > 0)

                frd = float(np.nanmax(tab[indx_ref]['mag']))
                if not np.isfinite(frd):
                    frd = None

                phot_summary = {'fileid_img': fileid_img, 'fileid_phot': fileid,
                                'fileid_template': tab.meta['template'], 'fileid_diffimg': tab.meta['diffimg'],
                                'targetid': targetid, 'obstime': tab.meta['obstime-bmjd'],
                                'photfilter': tab.meta['photfilter'], 'mag_raw': float(tab[indx_raw]['mag']),
                                'mag_raw_error': float(tab[indx_raw]['mag_error']),
                                'mag_sub': None if not any(indx_sub) else float(tab[indx_sub]['mag']),
                                'mag_sub_error': None if not any(indx_sub) else float(tab[indx_sub]['mag_error']),
                                'zeropoint': float(tab.meta['zp']), 'zeropoint_error': float(tab.meta['zp_error']),
                                'zeropoint_diff': float(tab.meta['zp_diff']), 'fwhm': float(tab.meta['fwhm'].value),
                                'seeing': float(tab.meta['seeing'].value), 'references_detected': int(np.sum(indx_ref)),
                                'used_for_epsf': int(np.sum(tab['used_for_epsf'])), 'faintest_reference_detected': frd,
                                'pipeline_version': tab.meta['version'], 'latest_version': new_version}

                db.cursor.execute("""INSERT INTO flows.photometry_details (
                    fileid_phot,
                    fileid_img,
                    fileid_template,
                    fileid_diffimg,
                    obstime_bmjd,
                    mag_raw,
                    mag_raw_error,
                    mag_sub,
                    mag_sub_error,
                    zeropoint,
                    zeropoint_error,
                    zeropoint_diff,
                    fwhm,
                    seeing,
                    references_detected,
                    used_for_epsf,
                    faintest_reference_detected,
                    pipeline_version
                ) VALUES (
                    %(fileid_phot)s,
                    %(fileid_img)s,
                    %(fileid_template)s,
                    %(fileid_diffimg)s,
                    %(obstime)s,
                    %(mag_raw)s,
                    %(mag_raw_error)s,
                    %(mag_sub)s,
                    %(mag_sub_error)s,
                    %(zeropoint)s,
                    %(zeropoint_error)s,
                    %(zeropoint_diff)s,
                    %(fwhm)s,
                    %(seeing)s,
                    %(references_detected)s,
                    %(used_for_epsf)s,
                    %(faintest_reference_detected)s,
                    %(pipeline_version)s
                );""", phot_summary)

                db.cursor.execute("SELECT * FROM flows.photometry_summary WHERE fileid_img=%s;", [fileid_img])
                if db.cursor.fetchone() is None:
                    db.cursor.execute("""INSERT INTO flows.photometry_summary (
                        fileid_phot,
                        fileid_img,
                        fileid_template,
                        fileid_diffimg,
                        targetid,
                        obstime,
                        photfilter,
                        mag_raw,
                        mag_raw_error,
                        mag_sub,
                        mag_sub_error,
                        pipeline_version,
                        latest_version
                    ) VALUES (
                        %(fileid_phot)s,
                        %(fileid_img)s,
                        %(fileid_template)s,
                        %(fileid_diffimg)s,
                        %(targetid)s,
                        %(obstime)s,
                        %(photfilter)s,
                        %(mag_raw)s,
                        %(mag_raw_error)s,
                        %(mag_sub)s,
                        %(mag_sub_error)s,
                        %(pipeline_version)s,
                        %(latest_version)s
                    );""", phot_summary)
                else:
                    db.cursor.execute("""UPDATE flows.photometry_summary SET
                        fileid_phot=%(fileid_phot)s,
                        targetid=%(targetid)s,
                        fileid_template=%(fileid_template)s,
                        fileid_diffimg=%(fileid_diffimg)s,
                        obstime=%(obstime)s,
                        photfilter=%(photfilter)s,
                        mag_raw=%(mag_raw)s,
                        mag_raw_error=%(mag_raw_error)s,
                        mag_sub=%(mag_sub)s,
                        mag_sub_error=%(mag_sub_error)s,
                        pipeline_version=%(pipeline_version)s,
                        latest_version=%(latest_version)s
                        WHERE fileid_img=%(fileid_img)s;""", phot_summary)

                # Update the photometry status to done:
                db.cursor.execute(
                    "UPDATE flows.photometry_status SET status='done' WHERE fileid=%(fileid_img)s AND status='ingest';",
                    phot_summary)

                if uploadlogid:
                    db.cursor.execute("UPDATE flows.uploadlog SET fileid=%s,status='ok' WHERE logid=%s;",
                                      [fileid, uploadlogid])

                db.conn.commit()

            except (KeyboardInterrupt, SystemExit):
                db.conn.rollback()
                if newpath is not None and os.path.isdir(os.path.dirname(newpath)):
                    shutil.rmtree(os.path.dirname(newpath))
                logger.warning("Stopped by user or system")
                return

            except Exception as e: # noqa: E722, pragma: no cover
                db.conn.rollback()
                if newpath is not None and os.path.isdir(os.path.dirname(newpath)):
                    shutil.rmtree(os.path.dirname(newpath))
                logger.exception("%s: Could not ingest file", fpath)
                if uploadlogid:
                    msg = str(e) if not hasattr(e, 'message') else str(e.message)
                    db.cursor.execute("UPDATE flows.uploadlog SET status=%s WHERE logid=%s;", [msg, uploadlogid])
                    db.conn.commit()

            else:
                # Set file and directory permissions:
                # TODO: Can this not be handled in a more elegant way?
                os.chmod(os.path.join(rootdir_archive, targetname), 0o2750)
                os.chmod(os.path.join(rootdir_archive, targetname, f'{fileid_img:05d}'), 0o2750)
                os.chmod(os.path.join(rootdir_archive, targetname, f'{fileid_img:05d}', f'v{new_version:02d}'), 0o2550)
                for f in os.listdir(os.path.dirname(newpath)):
                    os.chmod(os.path.join(os.path.dirname(newpath), f), 0o0440)

                logger.info("DELETE THE ORIGINAL FILE")
                if os.path.isfile(fpath):
                    os.remove(fpath)
                if uploadlogid:
                    db.cursor.execute("UPDATE flows.uploadlog SET uploadpath=NULL WHERE logid=%s;", [uploadlogid])
                    db.conn.commit()


# --------------------------------------------------------------------------------------------------
def cleanup_inbox():
    """
    Cleanup of inbox directory
    """
    rootdir_inbox = '/data/flows/archive/inbox'

    # Just a simple check to begin with:
    if not os.path.isdir(rootdir_inbox):
        raise FileNotFoundError("INBOX could not be found.")

    # Remove empty directories:
    for inputtype in ('science', 'templates', 'subtracted', 'photometry', 'replace'):
        for dpath in glob.iglob(os.path.join(rootdir_inbox, '*', inputtype)):
            if not os.listdir(dpath):
                os.rmdir(dpath)

    for dpath in glob.iglob(os.path.join(rootdir_inbox, '*')):
        if os.path.isdir(dpath) and not os.listdir(dpath):
            os.rmdir(dpath)

    # Delete left-over files in the database tables, that have been removed from disk:
    with AADC_DB() as db:
        db.cursor.execute("SELECT logid,uploadpath FROM flows.uploadlog WHERE uploadpath IS NOT NULL AND fileid IS NULL AND status!='File deleted';")
        for row in db.cursor.fetchall():
            if not os.path.isfile(os.path.join(rootdir_inbox, row['uploadpath'])):
                print(f"MARK AS DELETED IN DATABASE (logid={row['logid']:d}): " + row['uploadpath'])
                db.cursor.execute("UPDATE flows.uploadlog SET status='File deleted' WHERE logid=%s;", [row['logid']])
        db.conn.commit()


# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    logging_level = logging.INFO

    # Setup logging:
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    console = logging.StreamHandler(sys.stdout)
    console.setFormatter(formatter)
    logger = logging.getLogger(__name__)
    if not logger.hasHandlers():
        logger.addHandler(console)
    logger.setLevel(logging_level)

    # Add a CounterFilter to the logger, which will count the number of log-records
    # being passed through the logger. Can be used to count the number of errors/warnings:
    counter = CounterFilter()
    logger.addFilter(counter)

    # Run the ingests and cleanup:
    ingest_from_inbox()
    ingest_photometry_from_inbox()
    cleanup_inbox()

    # Last time we warned:
    report_errors = True
    try:
        with AADC_DB() as db:
            db.cursor.execute("SELECT * FROM flows.ingest_reporting WHERE last_ingest_report > TIMEZONE('utc', NOW()) - '1 day'::interval LIMIT 1;")
            last_ingest_report = db.cursor.fetchone()
            if last_ingest_report is None:
                db.cursor.execute("TRUNCATE flows.ingest_reporting;")
                db.cursor.execute("INSERT INTO flows.ingest_reporting (last_ingest_report) VALUES (TIMEZONE('utc', NOW()));")
                db.conn.commit()
                report_errors = True
            else:
                report_errors = False
    except:
        db.conn.rollback()
        report_errors = True

    # Check the number of errors or warnings issued, and convert these to a return-code:
    if report_errors:
        logcounts = counter.counter
        if logcounts.get('ERROR', 0) > 0 or logcounts.get('CRITICAL', 0) > 0:
            sys.exit(4)
        elif logcounts.get('WARNING', 0) > 0:
            sys.exit(3)
    sys.exit(0)
