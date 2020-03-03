#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

.. codeauthor:: Rasmus Handberg <rasmush@phys.au.dk>
"""

import logging
import os.path
import subprocess
import shlex
from .config import load_config
from .aadc_db import AADC_DB

#--------------------------------------------------------------------------------------------------
class CasjobsException(Exception):
	pass

#--------------------------------------------------------------------------------------------------
def floatval(value):
	return None if value == '' or value == '0' else float(value)

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
	except:
		os.remove(casjobs_config)

#--------------------------------------------------------------------------------------------------
def query_casjobs_refcat2(ra, dec, radius=24.0/60.0):
	"""
	Uses the CasJobs program to do a cone-search around the position.

	Parameters:
		ra (float):
		dec (float):
		radius (float, optional):

	Returns:
		list: List of dicts with the REFCAT2 information.

	.. codeauthor:: Rasmus Handberg <rasmush@phys.au.dk>
	"""

	logger = logging.getLogger(__name__)

	sql = "SELECT r.* FROM fGetNearbyObjEq({ra:f}, {dec:f}, {radius:f}) AS n INNER JOIN HLSP_ATLAS_REFCAT2.refcat2 AS r ON n.objid=r.objid ORDER BY n.distance;".format(
		ra=ra,
		dec=dec,
		radius=radius
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
		raise CasjobsException("ERROR detected in CasJobs: " + error_msg)

	if not results:
		raise CasjobsException("Could not find anything on CasJobs")

	return results

#--------------------------------------------------------------------------------------------------
def download_catalog(target=None):

	logger = logging.getLogger(__name__)

	with AADC_DB() as db:

		# Get the information about the target from the database:
		if target is None:
			db.cursor.execute("SELECT targetid,ra,decl FROM flows.targets WHERE catalog_downloaded=FALSE;")
		else:
			db.cursor.execute("SELECT targetid,ra,decl FROM flows.targets WHERE targetid=%s;", (target,))

		for row in db.cursor.fetchall():
			# The unique identifier of the target:
			targetid = int(row['targetid'])

			# Query the REFCAT2 catalog using CasJobs around the target position:
			results = query_casjobs_refcat2(row['ra'], row['decl'])

			# Insert the catalog into the local database:
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
				g_mag,
				r_mag,
				i_mag,
				z_mag,
				"J_mag",
				"H_mag",
				"K_mag")
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
				%(g_mag)s,
				%(r_mag)s,
				%(i_mag)s,
				%(z_mag)s,
				%(J_mag)s,
				%(H_mag)s,
				%(K_mag)s)
			ON CONFLICT DO NOTHING;""", results)
			logger.info("%d catalog entries inserted.", db.cursor.rowcount)

			# Mark the target that the catalog has been downloaded:
			db.cursor.execute("UPDATE flows.targets SET catalog_downloaded=TRUE WHERE targetid=%s;", (targetid,))
			db.conn.commit()