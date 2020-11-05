#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

.. codeauthor:: Rasmus Handberg <rasmush@phys.au.dk>
"""

import logging
import os.path
import subprocess
import shlex
import requests
import numpy as np
from astropy.coordinates import SkyCoord, Angle
from astropy import units as u
from .config import load_config
from .aadc_db import AADC_DB
from .ztf import query_ztf_id

#--------------------------------------------------------------------------------------------------
class CasjobsException(Exception):
	pass

#--------------------------------------------------------------------------------------------------
class CasjobsMemoryError(Exception):
	pass

#--------------------------------------------------------------------------------------------------
def floatval(value):
	return None if value == '' or value == 'NA' or value == '0' else float(value)

#--------------------------------------------------------------------------------------------------
def intval(value):
	return None if value == '' else int(value)

#--------------------------------------------------------------------------------------------------
def configure_casjobs(overwrite=False):
	"""
	Set up CasJobs if needed.

	Parameters:
		overwrite (bool, optional): Overwrite existing configuration. Default (False) is to not
			overwrite existing configuration.

	.. codeauthor:: Rasmus Handberg <rasmush@phys.au.dk>
	"""

	__dir__ = os.path.dirname(os.path.realpath(__file__))
	casjobs_config = os.path.join(__dir__, 'casjobs', 'CasJobs.config')
	if os.path.isfile(casjobs_config) and not overwrite:
		return

	config = load_config()
	wsid = config.get('casjobs', 'wsid', fallback=None)
	passwd = config.get('casjobs', 'password', fallback=None)
	if wsid is None or passwd is None:
		raise CasjobsException("CasJobs WSID and PASSWORD not in config.ini")

	try:
		with open(casjobs_config, 'w') as fid:
			fid.write("wsid={0:s}\n".format(wsid))
			fid.write("password={0:s}\n".format(passwd))
			fid.write("default_target=HLSP_ATLAS_REFCAT2\n")
			fid.write("default_queue=1\n")
			fid.write("default_days=1\n")
			fid.write("verbose=false\n")
			fid.write("debug=false\n")
			fid.write("jobs_location=http://mastweb.stsci.edu/gcasjobs/services/jobs.asmx\n")
	except: # noqa: E722, pragma: no cover
		if os.path.isfile(casjobs_config):
			os.remove(casjobs_config)

#--------------------------------------------------------------------------------------------------
def query_casjobs_refcat2(coo_centre, radius=24*u.arcmin):
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

	logger = logging.getLogger(__name__)
	if isinstance(radius, (float, int)):
		radius *= u.deg

	try:
		results = _query_casjobs_refcat2(coo_centre, radius=radius)
	except CasjobsMemoryError:
		logger.debug("CasJobs failed with memory error. Trying to use smaller radii.")
		results = _query_casjobs_refcat2_divide_and_conquer(coo_centre, radius=radius)

	# Remove duplicate entries:
	_, indx = np.unique([res['starid'] for res in results], return_index=True)
	results = [results[k] for k in indx]

	# Trim away anything outside radius:
	ra = [res['ra'] for res in results]
	decl = [res['decl'] for res in results]
	coords = SkyCoord(ra=ra, dec=decl, unit='deg', frame='icrs')
	sep = coords.separation(coo_centre)
	results = [res for k,res in enumerate(results) if sep[k] <= radius]

	logger.debug("Found %d unique results", len(results))
	return results

#--------------------------------------------------------------------------------------------------
def _query_casjobs_refcat2_divide_and_conquer(coo_centre, radius):
	logger = logging.getLogger(__name__)

	# Just put in a stop criterion to avoid infinite recursion:
	if radius < 0.04*u.deg:
		raise Exception("Too many subdivides")

	# Search central cone:
	try:
		results = _query_casjobs_refcat2(coo_centre, radius=0.5*radius)
	except CasjobsMemoryError:
		logger.debug("CasJobs failed with memory error. Trying to use smaller radii.")
		results = _query_casjobs_refcat2_divide_and_conquer(coo_centre, radius=0.5*radius)

	# Search six cones around central cone:
	for n in range(6):
		# FIXME: The 0.8 here is kind of a guess. There should be an analytic solution
		new = SkyCoord(
			ra=coo_centre.ra.deg + 0.8 * Angle(radius).deg * np.cos(n*60*np.pi/180),
			dec=coo_centre.dec.deg + 0.8 * Angle(radius).deg * np.sin(n*60*np.pi/180),
			unit='deg', frame='icrs')

		try:
			results += _query_casjobs_refcat2(new, radius=0.5*radius)
		except CasjobsMemoryError:
			logger.debug("CasJobs failed with memory error. Trying to use smaller radii.")
			results += _query_casjobs_refcat2_divide_and_conquer(new, radius=0.5*radius)

	return results

#--------------------------------------------------------------------------------------------------
def _query_casjobs_refcat2(coo_centre, radius=24*u.arcmin):
	"""
	Uses the CasJobs program to do a cone-search around the position.

	Parameters:
		coo_centre (:class:`astropy.coordinates.SkyCoord`): Coordinates of centre of search cone.
		radius (Angle, optional): Search radius. Default is 24 arcmin.

	Returns:
		list: List of dicts with the REFCAT2 information.

	.. codeauthor:: Rasmus Handberg <rasmush@phys.au.dk>
	"""

	logger = logging.getLogger(__name__)
	if isinstance(radius, (float, int)):
		radius *= u.deg

	sql = "SELECT r.* FROM fGetNearbyObjEq({ra:f}, {dec:f}, {radius:f}) AS n INNER JOIN HLSP_ATLAS_REFCAT2.refcat2 AS r ON n.objid=r.objid ORDER BY n.distance;".format(
		ra=coo_centre.ra.deg,
		dec=coo_centre.dec.deg,
		radius=Angle(radius).deg
	)
	logger.debug(sql)

	# Make sure that CasJobs have been configured:
	configure_casjobs()

	# The command to run the casjobs script:
	# BEWARE: This may change in the future without warning - it has before!
	cmd = 'java -jar casjobs.jar execute "{0:s}"'.format(sql)

	# Execute the command:
	cmd = shlex.split(cmd)
	directory = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'casjobs')
	proc = subprocess.Popen(cmd, cwd=directory, stdout=subprocess.PIPE, universal_newlines=True)
	stdout, stderr = proc.communicate()
	output = stdout.split("\n")

	# build list of all kois from output from the CasJobs-script:
	error_thrown = False
	results = []
	for line in output:
		line = line.strip()
		if line == '': continue
		if 'ERROR' in line:
			logger.error(line)
			error_thrown = True
			break

		row = line.split(',')
		if len(row) == 45 and row[0] != '[objid]:Integer':
			results.append({
				'starid': int(row[0]),
				'ra': floatval(row[1]),
				'decl': floatval(row[2]),
				'pm_ra': floatval(row[5]),
				'pm_dec': floatval(row[7]),
				'gaia_mag': floatval(row[9]),
				'gaia_bp_mag': floatval(row[11]),
				'gaia_rp_mag': floatval(row[13]),
				'gaia_variability': intval(row[17]),
				'g_mag': floatval(row[22]),
				'r_mag': floatval(row[26]),
				'i_mag': floatval(row[30]),
				'z_mag': floatval(row[34]),
				'J_mag': floatval(row[39]),
				'H_mag': floatval(row[41]),
				'K_mag': floatval(row[43]),
			})

	if error_thrown:
		error_msg = ''
		for line in output:
			if len(line.strip()) > 0:
				error_msg += line.strip() + "\n"

		logger.debug("Error Msg: %s", error_msg)
		if 'query results exceed memory limit' in error_msg.lower():
			raise CasjobsMemoryError("Query results exceed memory limit")
		else:
			raise CasjobsException("ERROR detected in CasJobs: " + error_msg)

	if not results:
		raise CasjobsException("Could not find anything on CasJobs")

	logger.debug("Found %d results", len(results))
	return results

#--------------------------------------------------------------------------------------------------
def query_apass(coo_centre, radius=24*u.arcmin):
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

	r = requests.post('https://www.aavso.org/cgi-bin/apass_dr10_download.pl',
		data={'ra': coo_centre.ra.deg, 'dec': coo_centre.dec.deg, 'radius': Angle(radius).deg, 'outtype': '1'})

	results = []

	lines = r.text.split("\n")
	#header = lines[0]

	for line in lines[1:]:
		if line.strip() == '': continue
		row = line.strip().split(',')

		results.append({
			'ra': floatval(row[0]),
			'decl': floatval(row[2]),
			'V_mag': floatval(row[4]),
			'B_mag': floatval(row[7]),
			'u_mag': floatval(row[10]),
			'g_mag': floatval(row[13]),
			'r_mag': floatval(row[16]),
			'i_mag': floatval(row[19]),
			'z_mag': floatval(row[22]),
			'Y_mag': floatval(row[25])
		})

	return results

#--------------------------------------------------------------------------------------------------
def query_all(coo_centre, radius=24*u.arcmin, dist_cutoff=2*u.arcsec):
	"""
	Query all catalogs, and return merged catalog.

	Parameters:
		coo_centre (:class:`astropy.coordinates.SkyCoord`): Coordinates of centre of search cone.
		radius (float): Search radius. Default 24 arcmin.
		dist_cutoff (float): Maximal distance between object is catalog matching. Default 2 arcsec.

	Returns:
		list: List of dicts with catalog stars.

	.. codeauthor:: Rasmus Handberg <rasmush@phys.au.dk>
	"""

	# Query the REFCAT2 catalog using CasJobs around the target position:
	results = query_casjobs_refcat2(coo_centre, radius=radius)

	# Query APASS around the target position:
	results_apass = query_apass(coo_centre, radius=radius)

	# Match the two catalogs using coordinates:
	# https://docs.astropy.org/en/stable/coordinates/matchsep.html#matching-catalogs
	ra = np.array([r['ra'] for r in results])
	decl = np.array([r['decl'] for r in results])
	refcat = SkyCoord(ra=ra, dec=decl, unit=u.deg, frame='icrs')

	ra_apass = np.array([r['ra'] for r in results_apass])
	decl_apass = np.array([r['decl'] for r in results_apass])
	apass = SkyCoord(ra=ra_apass, dec=decl_apass, unit=u.deg, frame='icrs')

	# Match the two catalogs:
	idx, d2d, _ = apass.match_to_catalog_sky(refcat)

	# Go through the matches and make sure they are valid:
	for k, i in enumerate(idx):
		# If APASS doesn't contain any new information anyway, skip it:
		if results_apass[k]['B_mag'] is None and results_apass[k]['V_mag'] is None \
			and results_apass[k]['u_mag'] is None:
			continue

		# Reject any match further away than the cutoff:
		if d2d[k] > dist_cutoff:
			continue

		# TODO: Use the overlapping magnitudes to make better match:
		#photdist = 0
		#for photfilt in ('g_mag', 'r_mag', 'i_mag', 'z_mag'):
		#	if results_apass[k][photfilt] and results[i][photfilt]:
		#		photdist += (results[i][photfilt] - results_apass[k][photfilt])**2
		#print( np.sqrt(photdist) )

		# Update the results "table" with the APASS filters:
		results[i].update({
			'V_mag': results_apass[k]['V_mag'],
			'B_mag': results_apass[k]['B_mag'],
			'u_mag': results_apass[k]['u_mag']
		})

	# Fill in empty fields where nothing was matched:
	for k in range(len(results)):
		if 'V_mag' not in results[k]:
			results[k].update({
				'B_mag': None,
				'V_mag': None,
				'u_mag': None
			})

	return results

#--------------------------------------------------------------------------------------------------
def download_catalog(target=None, radius=24*u.arcmin, dist_cutoff=2*u.arcsec):

	logger = logging.getLogger(__name__)

	with AADC_DB() as db:

		# Get the information about the target from the database:
		if target is not None and isinstance(target, int):
			db.cursor.execute("SELECT targetid,target_name,ra,decl FROM flows.targets WHERE targetid=%s;", [target])
		elif target is not None:
			db.cursor.execute("SELECT targetid,target_name,ra,decl FROM flows.targets WHERE target_name=%s;", [target])
		else:
			db.cursor.execute("SELECT targetid,target_name,ra,decl FROM flows.targets WHERE catalog_downloaded=FALSE;")

		for row in db.cursor.fetchall():
			# The unique identifier of the target:
			targetid = int(row['targetid'])
			target_name = row['target_name']

			# Coordinate of the target, which is the centre of the search cone:
			coo_centre = SkyCoord(ra=row['ra'], dec=row['decl'], unit=u.deg, frame='icrs')

			results = query_all(coo_centre, radius=radius, dist_cutoff=dist_cutoff)

			# Query for a ZTF identifier for this target:
			ztf_id = query_ztf_id(coo_centre, radius=radius)

			# Insert the catalog into the local database:
			try:
				#db.cursor.execute("TRUNCATE flows.refcat2;")
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
				ON CONFLICT DO NOTHING;""", results)
				logger.info("%d catalog entries inserted for %s.", db.cursor.rowcount, target_name)

				# Mark the target that the catalog has been downloaded:
				db.cursor.execute("UPDATE flows.targets SET catalog_downloaded=TRUE,ztf_id=%s WHERE targetid=%s;", (ztf_id, targetid))
				db.conn.commit()
			except:
				db.conn.rollback()
				raise
