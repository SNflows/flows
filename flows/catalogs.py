#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

.. codeauthor:: Rasmus Handberg <rasmush@phys.au.dk>
"""
import configparser
import logging
import os
import os.path
from pickle import NONE
import shlex
import subprocess
import warnings
from io import BytesIO

import numpy as np
import requests
from astropy import units as u
from astropy.coordinates import Angle, SkyCoord
from astropy.table import MaskedColumn, Table, unique
from astropy.time import Time
from astroquery import sdss
from astroquery.simbad import Simbad
from bottleneck import anynan
from tendrils.utils import load_config, query_ztf_id

import mastcasjobs

from .aadc_db import AADC_DB

logger = logging.getLogger(__name__)

# --------------------------------------------------------------------------------------------------
class CasjobsError(RuntimeError):
    pass


# --------------------------------------------------------------------------------------------------
class CasjobsMemoryError(RuntimeError):
    pass


# --------------------------------------------------------------------------------------------------
def floatval(value):
    return None if value == '' or value == 'NA' or value == '0' else float(value)


# --------------------------------------------------------------------------------------------------
def intval(value):
    return None if value == '' else int(value)


# --------------------------------------------------------------------------------------------------
def casjobs_configured(config: configparser.ConfigParser) -> bool:
    """
    Check if casjobs credentials are available.
    """
    wsid = config.get('casjobs', 'wsid', fallback=os.environ.get("CASJOBS_WSID", None))
    passwd = config.get('casjobs', 'password', fallback=os.environ.get("CASJOBS_PASSWORD", None))
    if wsid is None or passwd is None:
        raise CasjobsError("CasJobs WSID and PASSWORD not in config.ini / 'CASJOBS_WSID' and 'CASJOBS_PASSWORD' not defined as environment variables")
    return True

# --------------------------------------------------------------------------------------------------
def query_casjobs_refcat2(coo_centre, radius=24 * u.arcmin):
    """
    Uses the CasJobs program to do a cone-search around the position.

    Will first attempt to do single large cone-search, and if that
    fails because of CasJobs memory limits, will sub-divide the cone
    into smaller queries.

    Parameters:
        coo_centre (:class:`astropy.coordinates.SkyCoord`): Coordinates of centre of search cone.
        radius (Angle, optional): Search radius. Default is 24 arcmin.

    Returns:
        list: List of dicts with the REFCAT2 information.

    .. codeauthor:: Rasmus Handberg <rasmush@phys.au.dk>
    """
    if isinstance(radius, (float, int)):
        radius *= u.deg

    try:
        results = _query_casjobs_refcat2(coo_centre, radius=radius)
    except CasjobsMemoryError:
        logger.debug("CasJobs failed with memory error. Trying to use smaller radii.")
        results = _query_casjobs_refcat2_divide_and_conquer(coo_centre, radius=radius)

    # Remove duplicate entries:
    results = unique(results, keys='starid')

    # Trim away anything outside radius:
    ra = [res['ra'] for res in results]
    decl = [res['decl'] for res in results]
    coords = SkyCoord(ra=ra, dec=decl, unit='deg', frame='icrs')
    sep = coords.separation(coo_centre)
    mask = sep <= radius

    logger.debug("Found %d unique results", np.count_nonzero(mask))
    return results[mask]


# --------------------------------------------------------------------------------------------------
def _query_casjobs_refcat2_divide_and_conquer(coo_centre, radius):

    # Just put in a stop criterion to avoid infinite recursion:
    if radius < 0.04 * u.deg:
        raise Exception("Too many subdivides")

    # Search central cone:
    try:
        results = _query_casjobs_refcat2(coo_centre, radius=0.5 * radius)
    except CasjobsMemoryError:
        logger.debug("CasJobs failed with memory error. Trying to use smaller radii.")
        results = _query_casjobs_refcat2_divide_and_conquer(coo_centre, radius=0.5 * radius)

    # Search six cones around central cone:
    for n in range(6):
        # FIXME: The 0.8 here is kind of a guess. There should be an analytic solution
        new = SkyCoord(ra=coo_centre.ra.deg + 0.8 * Angle(radius).deg * np.cos(n * 60 * np.pi / 180),
                       dec=coo_centre.dec.deg + 0.8 * Angle(radius).deg * np.sin(n * 60 * np.pi / 180), unit='deg',
                       frame='icrs')

        try:
            results += _query_casjobs_refcat2(new, radius=0.5 * radius)
        except CasjobsMemoryError:
            logger.debug("CasJobs failed with memory error. Trying to use smaller radii.")
            results += _query_casjobs_refcat2_divide_and_conquer(new, radius=0.5 * radius)

    return results


# --------------------------------------------------------------------------------------------------
def _query_casjobs_refcat2(coo_centre, radius=24 * u.arcmin):
    """
    Uses the CasJobs program to do a cone-search around the position.

    Parameters:
        coo_centre (:class:`astropy.coordinates.SkyCoord`): Coordinates of centre of search cone.
        radius (Angle, optional): Search radius. Default is 24 arcmin.

    Returns:
        list: List of dicts with the REFCAT2 information.

    .. codeauthor:: Rasmus Handberg <rasmush@phys.au.dk>
    """

    if isinstance(radius, (float, int)):
        radius *= u.deg

    # prepare refcat2 query to find all objects near position
    casjobs_context = "HLSP_ATLAS_REFCAT2"
    casjobs_personal_table = "flows"
    sql = ("""SELECT r.* INTO {table:s} FROM fGetNearbyObjEq({ra:f}, {dec:f}, {radius:f}) AS n
    INNER JOIN {context:s}.refcat2 AS r ON n.objid=r.objid ORDER BY n.distance;"""
    .format(table=casjobs_personal_table, context=casjobs_context,
            ra=coo_centre.ra.deg, dec=coo_centre.dec.deg, radius=Angle(radius).deg))

    config = load_config()
    casjobs_configured(config)
    jobs = mastcasjobs.MastCasJobs(userid=config.get('casjobs', 'wsid'),
                                   password=config.get('casjobs', 'password'),
                                   context=casjobs_context)

    # limited storage space at remote casjobs, drop table, if it exists, before making a new one
    for table in jobs.list_tables():
        if table == casjobs_personal_table:
            jobs.drop_table(casjobs_personal_table)  # more graceful than drop_table_if_exists()

    # submit job, wait for it to finish and retrieve the results
    job_id = jobs.submit(sql, context=casjobs_context, task_name=casjobs_personal_table)
    logger.debug("Submitted query '%s' to MAST CasJobs context '%s' with task name '%s'",
                 sql, casjobs_context, casjobs_personal_table)
    code, status = jobs.monitor(job_id, timeout=5)
    logger.debug("MAST CasJobs query exited with status %d: %s", code, status)
    results = jobs.get_table(casjobs_personal_table)

    # Contents of HLSP_ATLAS_REFCAT2.refcat2 data table - i.e. the columns of 'results' just
    # above; from https://galex.stsci.edu/casjobs/MyDB.aspx -> HLSP_ATLAS_REFCAT2.refcat2
    # -----------------------------------------------------------------------------
    # Name      Unit           Data Type  Size  Default Value  Description
    # objid     dimensionless  bigint     8                    Unique object identifier.
    # RA        degrees        float      8                    RA from Gaia DR2, J2000, epoch 2015.5
    # Dec       degrees        float      8                    Dec from Gaia DR2, J2000, epoch 2015.5
    # plx       mas            real       4                    Parallax from Gaia DR2
    # dplx      mas            real       4                    Parallax uncertainty
    # pmra      mas/yr         real       4                    Proper motion in RA from Gaia DR2
    # dpmra     mas/yr         real       4                    Proper motion uncertainty in RA
    # pmdec     mas/yr         real       4                    Proper motion in Dec from Gaia DR2
    # dpmdec    mas/yr         real       4                    Proper motion uncertainty in Dec
    # Gaia      mag            real       4                    Gaia DR2 G magnitude
    # dGaia     mag            real       4                    Gaia DR2 G magnitude uncertainty
    # BP        mag            real       4                    Gaia G_BP magnitude
    # dBP       mag            real       4                    Gaia G_BP magnitude uncertainty
    # RP        mag            real       4                    Gaia G_RP magnitude
    # dRP       mag            real       4                    Gaia G_RP magnitude uncertainty
    # Teff      [K]            int        4                    Gaia stellar effective temperature
    # AGaia     mag            real       4                    Gaia estimate of G-band extinction for this star
    # dupvar    dimensionless  int        4                    Gaia flags coded as CONSTANT (0), VARIABLE (1), or NOT_AVAILABLE (2) + 4*DUPLICATE
    # Ag        mag            real       4                    SFD estimate of total column g-band extinction
    # rp1       arc seconds    real       4                    Radius where cumulative G flux exceeds 0.1x this star
    # r1        arc seconds    real       4                    Radius where cumulative G flux exceeds 1x this star
    # r10       arc seconds    real       4                    Radius where cumulative G flux exceeds 10x this star
    # g         mag            real       4                    Pan-STARRS g_P1 magnitude
    # dg        mag            real       4                    Pan-STARRS g_P1 magnitude uncertainty
    # gchi      dimensionless  real       4                    chi^2/DOF for contributors to g
    # gcontrib  dimensionless  int        4                    Bitmap of contributing catalogs to g
    # r         mag            real       4                    Pan-STARRS r_P1 magnitude
    # dr        mag            real       4                    Pan-STARRS r_P1 magnitude uncertainty
    # rchi      dimensionless  real       4                    chi^2/DOF for contributors to r
    # rcontrib  dimensionless  int        4                    Bitmap of contributing catalogs to r
    # i         mag            real       4                    Pan-STARRS i_P1 magnitude
    # di        mag            real       4                    Pan-STARRS i_P1 magnitude uncertainty
    # ichi      dimensionless  real       4                    chi^2/DOF for contributors to i
    # icontrib  dimensionless  int        4                    Bitmap of contributing catalogs to i
    # z         mag            real       4                    Pan-STARRS z_P1 magnitude
    # dz        mag            real       4                    Pan-STARRS z_P1 magnitude uncertainty
    # zchi      dimensionless  real       4                    chi^2/DOF for contributors to z
    # zcontrib  dimensionless  int        4                    Bitmap of contributing catalogs to z
    # nstat     dimensionless  int        4                    Count of griz deweighted outliers
    # J         mag            real       4                    2MASS J magnitude
    # dJ        mag            real       4                    2MASS J magnitude uncertainty
    # H         mag            real       4                    2MASS H magnitude
    # dH        mag            real       4                    2MASS H magnitude uncertainty
    # K         mag            real       4                    2MASS K magnitude
    # dK        mag            real       4                    2MASS K magnitude uncertainty

    # remove unused columns from results
    results.remove_columns(['plx', 'dplx', 'dpmra', 'dpmdec', 'dGaia', 'dBP', 'dRP', 'Teff', 'AGaia', 'Ag',
                            'rp1', 'r1', 'r10', 'dg', 'gchi', 'gcontrib', 'dr', 'rchi', 'rcontrib',
                            'di', 'ichi', 'icontrib', 'dz', 'zchi', 'zcontrib', 'nstat', 'dJ', 'dH', 'dK'])

    # rename columns to match older java-based implementation which had different naming conventions
    # TODO: refactor dependent functions to expect the default namings and get rid of these renamings
    #       to start: comment out the renamings below and see where the flows pipeline breaks, then fix it
    results = Table(results,
                    names=('starid', 'ra', 'decl', 'pm_ra', 'pm_dec', 'gaia_mag', 'gaia_bp_mag', 'gaia_rp_mag',
                           'gaia_variability', 'g_mag', 'r_mag', 'i_mag', 'z_mag', 'J_mag', 'H_mag', 'K_mag'))

    if not results:
        raise CasjobsError("Could not find anything on CasJobs")

    logger.debug("Found %d results", len(results))
    return results


# --------------------------------------------------------------------------------------------------
def query_apass(coo_centre, radius=24 * u.arcmin):
    """
    Queries APASS catalog using cone-search around the position.

    Parameters:
        coo_centre (:class:`astropy.coordinates.SkyCoord`): Coordinates of centre of search cone.
        radius (float, optional):

    Returns:
        list: List of dicts with the APASS information.

    .. codeauthor:: Rasmus Handberg <rasmush@phys.au.dk>
    """

    # https://vizier.u-strasbg.fr/viz-bin/VizieR-3?-source=II/336

    if isinstance(radius, (float, int)):
        radius *= u.deg

    data = {'ra': coo_centre.icrs.ra.deg, 'dec': coo_centre.icrs.dec.deg, 'radius': Angle(radius).deg, 'outtype': '1'}

    res = requests.post('https://www.aavso.org/cgi-bin/apass_dr10_download.pl', data=data)
    res.raise_for_status()

    results = []

    lines = res.text.split("\n")
    # header = lines[0]

    for line in lines[1:]:
        if line.strip() == '': continue
        row = line.strip().split(',')

        results.append(
            {'ra': floatval(row[0]), 'decl': floatval(row[2]), 'V_mag': floatval(row[4]), 'B_mag': floatval(row[7]),
             'u_mag': floatval(row[10]), 'g_mag': floatval(row[13]), 'r_mag': floatval(row[16]),
             'i_mag': floatval(row[19]), 'z_mag': floatval(row[22]), 'Y_mag': floatval(row[25])})

    return results


# --------------------------------------------------------------------------------------------------
def query_sdss(coo_centre, radius=24 * u.arcmin, dr=16, clean=True):
    """
    Queries SDSS catalog using cone-search around the position.

    Parameters:
        coo_centre (:class:`astropy.coordinates.SkyCoord`): Coordinates of centre of search cone.
        radius (float, optional):
        dr (int, optional): SDSS Data Release to query. Default=16.
        clean (bool, optional): Clean results for stars only and no other problems.

    Returns:
        tuple:
        - :class:`astropy.table.Table`: Table with SDSS information.
        - :class:`astropy.coordinates.SkyCoord`: Sky coordinates for SDSS objects.

    .. codeauthor:: Emir Karamehmetoglu <emir.k@phys.au.dk>
    .. codeauthor:: Rasmus Handberg <rasmush@phys.au.dk>
    """

    if isinstance(radius, (float, int)):
        radius *= u.deg

    # from astroquery version 0.4.7 onwards, radius is limited to 3 arcmins - use this value if necessary
    radius = radius if radius < sdss.SDSS.MAX_CROSSID_RADIUS else sdss.SDSS.MAX_CROSSID_RADIUS

    sdss.conf.skyserver_baseurl = sdss.conf.skyserver_baseurl.replace("http://","https://")
    AT_sdss = sdss.SDSS.query_region(coo_centre, photoobj_fields=['type', 'clean', 'ra', 'dec', 'psfMag_u'], data_release=dr,
                                timeout=600, radius=radius)

    if AT_sdss is None:
        return None, None

    if clean:
        # Clean SDSS following https://www.sdss.org/dr12/algorithms/photo_flags_recommend/
        # 6 == star, clean means remove interp, edge, suspicious defects, deblending problems, duplicates.
        AT_sdss = AT_sdss[(AT_sdss['type'] == 6) & (AT_sdss['clean'] == 1)]

    # Remove these columns since they are no longer needed:
    AT_sdss.remove_columns(['type', 'clean'])

    if len(AT_sdss) == 0:
        return None, None

    # Create SkyCoord object with the coordinates:
    sdss_coord = SkyCoord(ra=AT_sdss['ra'], dec=AT_sdss['dec'], unit=u.deg, frame='icrs')

    return AT_sdss, sdss_coord


# --------------------------------------------------------------------------------------------------
def query_simbad(coo_centre, radius=24 * u.arcmin):
    """
    Query SIMBAD using cone-search around the position using astroquery.

    Parameters:
        coo_centre (:class:`astropy.coordinates.SkyCoord`): Coordinates of centre of search cone.
        radius (float, optional):

    Returns:
        list: Astropy Table with SIMBAD information.

    .. codeauthor:: Rasmus Handberg <rasmush@phys.au.dk>
    """

    s = Simbad()
    s.ROW_LIMIT = 0
    s.add_votable_fields('ra(d;A;ICRS;J2000)', 'dec(d;D;ICRS;2000)', 'pmra', 'pmdec')
    s.add_votable_fields('otype')
    s.add_votable_fields('flux(B)', 'flux(V)', 'flux(R)', 'flux(I)', 'flux(J)', 'flux(H)', 'flux(K)')
    s.add_votable_fields('flux(u)', 'flux(g)', 'flux(r)', 'flux(i)', 'flux(z)')

    rad = Angle(radius).arcmin
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore', category=UserWarning)
        results = s.query_criteria(
            f'region(circle, icrs, {coo_centre.icrs.ra.deg:.10f} {coo_centre.icrs.dec.deg:+.10f}, {rad}m)',
            otypes='Star')

    if not results:
        return None, None

    # Rename columns:
    results.rename_column('MAIN_ID', 'main_id')
    results.rename_column('RA_d_A_ICRS_J2000', 'ra')
    results.rename_column('DEC_d_D_ICRS_2000', 'dec')
    results.rename_column('PMRA', 'pmra')
    results.rename_column('PMDEC', 'pmdec')
    results.rename_column('FLUX_B', 'B_mag')
    results.rename_column('FLUX_V', 'V_mag')
    results.rename_column('FLUX_R', 'R_mag')
    results.rename_column('FLUX_I', 'I_mag')
    results.rename_column('FLUX_J', 'J_mag')
    results.rename_column('FLUX_H', 'H_mag')
    results.rename_column('FLUX_K', 'K_mag')
    results.rename_column('FLUX_u', 'u_mag')
    results.rename_column('FLUX_g', 'g_mag')
    results.rename_column('FLUX_r', 'r_mag')
    results.rename_column('FLUX_i', 'i_mag')
    results.rename_column('FLUX_z', 'z_mag')
    results.rename_column('OTYPE', 'otype')
    results.remove_column('SCRIPT_NUMBER_ID')
    results.sort(['V_mag', 'B_mag', 'H_mag'])

    # Filter out object types which shouldn'r really be in there anyway:
    indx = (results['otype'] == 'Galaxy') | (results['otype'] == 'LINER') | (results['otype'] == 'SN')
    results = results[~indx]

    if len(results) == 0:
        return None, None

    # Build sky coordinates object:
    simbad = SkyCoord(ra=results['ra'], dec=results['dec'], pm_ra_cosdec=results['pmra'], pm_dec=results['pmdec'],
                      frame='icrs', obstime='J2000')

    return results, simbad


# --------------------------------------------------------------------------------------------------
def query_skymapper(coo_centre, radius=24 * u.arcmin):
    """
    Queries SkyMapper catalog using cone-search around the position.

    Parameters:
        coo_centre (:class:`astropy.coordinates.SkyCoord`): Coordinates of centre of search cone.
        radius (float, optional):

    Returns:
        tuple:
        - :class:`astropy.table.Table`: Astropy Table with SkyMapper information.
        - :class:`astropy.coordinates.SkyCoord`:

    .. codeauthor:: Rasmus Handberg <rasmush@phys.au.dk>
    """

    if isinstance(radius, (float, int)):
        radius *= u.deg

    # Query the SkyMapper cone-search API:
    params = {'RA': coo_centre.icrs.ra.deg, 'DEC': coo_centre.icrs.dec.deg, 'SR': Angle(radius).deg,
              'CATALOG': 'dr2.master', 'VERB': 1, 'RESPONSEFORMAT': 'VOTABLE'}
    res = requests.get('http://skymapper.anu.edu.au/sm-cone/public/query', params=params)
    res.raise_for_status()

    # For some reason the VOTable parser needs a file-like object:
    fid = BytesIO(bytes(res.text, 'utf8'))
    results = Table.read(fid, format='votable')

    if len(results) == 0:
        return None, None

    # Clean the results:
    # http://skymapper.anu.edu.au/data-release/dr2/#Access
    indx = (results['flags'] == 0) & (results['nimaflags'] == 0) & (results['ngood'] > 1)
    results = results[indx]
    if len(results) == 0:
        return None, None

    # Create SkyCoord object containing SkyMapper objects with their observation time:
    skymapper = SkyCoord(ra=results['raj2000'], dec=results['dej2000'],
                         obstime=Time(results['mean_epoch'], format='mjd', scale='utc'), frame='icrs')

    return results, skymapper


# --------------------------------------------------------------------------------------------------

def query_local_refcat2(coo_centre, radius=24 * u.arcmin):
    """
    Uses refcat.c from https://archive.stsci.edu/hlsp/atlas-refcat2 to do a cone-search around the position.
    NOTE: In going from mast casjobs to local refcat.c, we have lost the unique object identifier, 'objid'
               ----> See _query_casjobs_refcat2
          We have to invent our own objids.
          They get assigned incrementally from the value "300,000,000,000,000,000" onwards, which is comfortably
          above any objids in our database on 2025-01-25.
          TODO: Something less hacky! Assign every entry in the local catalog a unique objid somehow...

    Parameters:
        coo_centre (:class:`astropy.coordinates.SkyCoord`): Coordinates of centre of search cone.
        radius (Angle, optional): Search radius. Default is 24 arcmin.

    Returns:
        list: List of dicts with the REFCAT2 information.

    .. codeauthor:: Erik Jensen
    """

    current_max_objid = 300_000_000_000_000_000
    with AADC_DB() as db:
        db.cursor.execute("SELECT MAX(starid) FROM flows.refcat2;")
        for row in db.cursor.fetchall():
            if row[0] > current_max_objid:
                current_max_objid = row[0]

    refcat_output = subprocess.run(["refcat",
                              str(coo_centre.ra.deg), str(coo_centre.dec.deg),
                              "-rad", str(Angle(radius).deg),
                              "-var", "ra,dec,pmra,pmdec,gaia,bp,rp,dupvar,g,r,i,z,j,h,k",
                              "-nohdr"],
                             encoding="utf-8", capture_output=True, check=True)

    # TODO: feel free to make this less awful!
    OBJID = RA = DEC = PMRA = PMDEC = GAIA = BP = RP = DUPVAR = G = R = I = Z = J = H = K = np.array([])
    for line in refcat_output.stdout.splitlines():
        OBJID = np.append(OBJID, current_max_objid)
        current_max_objid += 1

        ra, dec, pmra, pmdec, gaia, bp, rp, dupvar, g, r, i, z, j, h, k = line.split()

        RA    = np.append(RA, ra)
        DEC   = np.append(DEC, dec)
        PMRA  = np.append(PMRA, pmra)
        PMDEC = np.append(PMDEC, pmdec)
        GAIA = np.append(GAIA, gaia)
        BP = np.append(BP, bp)
        RP = np.append(RP, rp)
        DUPVAR = np.append(DUPVAR, dupvar)
        G = np.append(G, g)
        R = np.append(R, r)
        I = np.append(I, i)
        Z = np.append(Z, z)
        J = np.append(J, j)
        H = np.append(H, h)
        K = np.append(K, k)

    # rename columns to match older java-based implementation which had different naming conventions
    # TODO: refactor dependent functions to expect the default namings and get rid of these renamings
    #       to start: comment out the renamings below and see where the flows pipeline breaks, then fix it
    results = Table({'starid': OBJID, 'ra': RA, 'decl': DEC, 'pm_ra': PMRA, 'pm_dec': PMDEC, 'gaia_mag': GAIA,
                'gaia_bp_mag': BP, 'gaia_rp_mag': RP, 'gaia_variability': DUPVAR, 'g_mag': G, 'r_mag': R, 'i_mag': I,
                'z_mag': Z, 'J_mag': J, 'H_mag': H, 'K_mag': K})

    if not results:
        raise CasjobsError("Something went wrong in local refcat.c binary")

    logger.debug("Found %d results", len(results))
    return results

def query_all(coo_centre, radius=24 * u.arcmin, dist_cutoff=2 * u.arcsec):
    """
    Query all catalogs (REFCAT2, APASS, SDSS and SkyMapper) and return merged catalog.

    Merging of catalogs are done using sky coordinates:
    https://docs.astropy.org/en/stable/coordinates/matchsep.html#matching-catalogs

    Parameters:
        coo_centre (:class:`astropy.coordinates.SkyCoord`): Coordinates of centre of search cone.
        radius (float): Search radius. Default 24 arcmin.
        dist_cutoff (float): Maximal distance between object is catalog matching. Default 2 arcsec.

    Returns:
        :class:`astropy.table.Table`: Table with catalog stars.

    TODO:
        - Use the overlapping magnitudes to make better matching.

    .. codeauthor:: Rasmus Handberg <rasmush@phys.au.dk>
    .. codeauthor:: Emir Karamehmetoglu <emir.k@phys.au.dk>
    """

    # Query the REFCAT2 catalog using CasJobs around the target position:
    results = query_local_refcat2(coo_centre, radius=radius)
    AT_results = Table(results)
    refcat = SkyCoord(ra=AT_results['ra'], dec=AT_results['decl'], unit=u.deg, frame='icrs')

    # REFCAT results table does not have uBV
    N = len(AT_results)
    d = np.full(N, np.nan)
    AT_results.add_column(MaskedColumn(name='B_mag', unit='mag', dtype='float32', fill_value=np.nan, data=d))
    AT_results.add_column(MaskedColumn(name='V_mag', unit='mag', dtype='float32', fill_value=np.nan, data=d))
    AT_results.add_column(MaskedColumn(name='u_mag', unit='mag', dtype='float32', fill_value=np.nan, data=d))

    # Query APASS around the target position:
    results_apass = query_apass(coo_centre, radius=radius)
    if results_apass:
        AT_apass = Table(results_apass)
        apass = SkyCoord(ra=AT_apass['ra'], dec=AT_apass['decl'], unit=u.deg, frame='icrs')

        # Match the two catalogs using coordinates:
        idx, d2d, _ = apass.match_to_catalog_sky(refcat)
        sep_constraint = (d2d <= dist_cutoff)  # Reject any match further away than the cutoff
        idx_apass = np.arange(len(idx), dtype='int')  # since idx maps apass to refcat

        # Update results table with APASS bands of interest
        AT_results['B_mag'][idx[sep_constraint]] = AT_apass[idx_apass[sep_constraint]]['B_mag']
        AT_results['V_mag'][idx[sep_constraint]] = AT_apass[idx_apass[sep_constraint]]['V_mag']
        AT_results['u_mag'][idx[sep_constraint]] = AT_apass[idx_apass[sep_constraint]]['u_mag']

    # Create SDSS around the target position:
    AT_sdss, sdss = query_sdss(coo_centre, radius=radius)
    if AT_sdss:
        # Match to dist_cutoff sky distance (angular) apart
        idx, d2d, _ = sdss.match_to_catalog_sky(refcat)
        sep_constraint = (d2d <= dist_cutoff)
        idx_sdss = np.arange(len(idx), dtype='int')  # since idx maps sdss to refcat

        # Overwrite APASS u-band with SDSS u-band:
        AT_results['u_mag'][idx[sep_constraint]] = AT_sdss[idx_sdss[sep_constraint]]['psfMag_u']

    # Query SkyMapper around the target position, only if there are missing u-band magnitudes:
    if anynan(AT_results['u_mag']):
        results_skymapper, skymapper = query_skymapper(coo_centre, radius=radius)
        if results_skymapper:
            idx, d2d, _ = skymapper.match_to_catalog_sky(refcat)
            sep_constraint = (d2d <= dist_cutoff)
            idx_skymapper = np.arange(len(idx), dtype='int')  # since idx maps skymapper to refcat

            newval = results_skymapper[idx_skymapper[sep_constraint]]['u_psf']
            oldval = AT_results['u_mag'][idx[sep_constraint]]
            indx = ~np.isfinite(oldval)
            if np.any(indx):
                AT_results['u_mag'][idx[sep_constraint]][indx] = newval[indx]

    return AT_results


# --------------------------------------------------------------------------------------------------
def convert_table_to_dict(tab):
    """
    Utility function for converting Astropy Table to list of dicts that the database
    likes as input.

    Parameters:
        tab (:class:`astropy.table.Table`): Astropy table coming from query_all.

    Returns:
        list: List of dicts where the column names are the keys. Datatypes are changed
            to things that the database understands (e.g. NaN -> None).

    .. codeauthor:: Rasmus Handberg <rasmush@phys.au.dk>
    """
    results = [dict(zip(tab.colnames, row)) for row in tab.filled()]
    for row in results:
        for key, val in row.items():
            if isinstance(val, (np.int64, np.int32)):
                row[key] = int(val)
            elif isinstance(val, (float, np.float32, np.float64)):
                if np.isfinite(val):
                    row[key] = float(val)
                else:
                    row[key] = None

    return results


# --------------------------------------------------------------------------------------------------
def download_catalog(target=None, radius=24 * u.arcmin, radius_ztf=3 * u.arcsec, dist_cutoff=2 * u.arcsec,
                     update_existing=False):
    """
    Download reference star catalogs and save to Flows database.

    Parameters:
        target (str or int): Target identifier to download catalog for.
        radius (Angle, optional): Radius around target to download catalogs.
        radius_ztf (Angle, optional): Radius around target to search for ZTF identifier.
        dist_cutoff (Angle, optional): Distance cutoff used for matching catalog positions.
        update_existing (bool, optional): Update existing catalog entries or skip them.

    .. codeauthor:: Rasmus Handberg <rasmush@phys.au.dk>
    """

    with AADC_DB() as db:

        # Get the information about the target from the database:
        if target is not None and isinstance(target, (int, float)):
            db.cursor.execute(
                "SELECT targetid,target_name,ra,decl,discovery_date FROM flows.targets WHERE targetid=%s;",
                [int(target)])
        elif target is not None:
            db.cursor.execute(
                "SELECT targetid,target_name,ra,decl,discovery_date FROM flows.targets WHERE target_name=%s;", [target])
        else:
            db.cursor.execute(
                "SELECT targetid,target_name,ra,decl,discovery_date FROM flows.targets WHERE catalog_downloaded=FALSE;")

        for row in db.cursor.fetchall():
            # The unique identifier of the target:
            targetid = int(row['targetid'])
            target_name = row['target_name']
            dd = row['discovery_date']
            if dd is not None:
                dd = Time(dd, format='datetime', scale='utc')

            # Coordinate of the target, which is the centre of the search cone:
            coo_centre = SkyCoord(ra=row['ra'], dec=row['decl'], unit=u.deg, frame='icrs')

            # Download combined catalog from all sources:
            tab = query_all(coo_centre, radius=radius, dist_cutoff=dist_cutoff)

            # Query for a ZTF identifier for this target:
            ztf_id = query_ztf_id(coo_centre, radius=radius_ztf, discovery_date=dd)

            # Because the database is picky with datatypes, we need to change things
            # before they are passed on to the database:
            results = convert_table_to_dict(tab)

            # Insert the catalog into the local database:
            if update_existing:
                on_conflict = """ON CONSTRAINT refcat2_pkey DO UPDATE SET
					ra=EXCLUDED.ra,
					decl=EXCLUDED.decl,
					pm_ra=EXCLUDED.pm_ra,
					pm_dec=EXCLUDED.pm_dec,
					gaia_mag=EXCLUDED.gaia_mag,
					gaia_bp_mag=EXCLUDED.gaia_bp_mag,
					gaia_rp_mag=EXCLUDED.gaia_rp_mag,
					gaia_variability=EXCLUDED.gaia_variability,
					u_mag=EXCLUDED.u_mag,
					g_mag=EXCLUDED.g_mag,
					r_mag=EXCLUDED.r_mag,
					i_mag=EXCLUDED.i_mag,
					z_mag=EXCLUDED.z_mag,
					"J_mag"=EXCLUDED."J_mag",
					"H_mag"=EXCLUDED."H_mag",
					"K_mag"=EXCLUDED."K_mag",
					"V_mag"=EXCLUDED."V_mag",
					"B_mag"=EXCLUDED."B_mag"
				WHERE refcat2.starid=EXCLUDED.starid"""
            else:
                on_conflict = 'DO NOTHING'

            # Avoid testing "ON CONFLICT" of postgres. Only test update/insert.
            if update_existing:
                try:
                    db.cursor.executemany("""INSERT INTO flows.refcat2 (
                        starid,
                        ra,
                        decl,
                        pm_ra,
                        pm_dec,
                        gaia_mag,
                        gaia_bp_mag,
                        gaia_rp_mag,
                        gaia_variability,
                        u_mag,
                        g_mag,
                        r_mag,
                        i_mag,
                        z_mag,
                        "J_mag",
                        "H_mag",
                        "K_mag",
                        "V_mag",
                        "B_mag")
                    VALUES (
                        %(starid)s,
                        %(ra)s,
                        %(decl)s,
                        %(pm_ra)s,
                        %(pm_dec)s,
                        %(gaia_mag)s,
                        %(gaia_bp_mag)s,
                        %(gaia_rp_mag)s,
                        %(gaia_variability)s,
                        %(u_mag)s,
                        %(g_mag)s,
                        %(r_mag)s,
                        %(i_mag)s,
                        %(z_mag)s,
                        %(J_mag)s,
                        %(H_mag)s,
                        %(K_mag)s,
                        %(V_mag)s,
                        %(B_mag)s)
                    ON CONFLICT """ + on_conflict + ";", results)
                    logger.info("%d catalog entries inserted for %s.", db.cursor.rowcount, target_name)

                    # Mark the target that the catalog has been downloaded:
                    db.cursor.execute("UPDATE flows.targets SET catalog_downloaded=TRUE,ztf_id=%s WHERE targetid=%s;",
                                    (ztf_id, targetid))
                    db.conn.commit()
                except:  # noqa: E722, pragma: no cover
                    db.conn.rollback()
                    raise
