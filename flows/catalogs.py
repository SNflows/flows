#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

.. codeauthor:: Rasmus Handberg <rasmush@phys.au.dk>
"""

import astropy.units as u
from astropy.table import Table
import os.path
import subprocess
import shlex
import logging
import requests
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

#--------------------------------------------------------------------------------------------------
def get_catalog(target, radius=None, output='table'):
	"""

	Parameters:
		target (int or str):
		radius (float, optional): Radius around target in degrees to return targets for.
		outout (str, optional): Desired output format. Choises are 'table', 'dict', 'json'.
			Default='table'.

	Returns:
		dict: Dictionary with three members:
			- 'target': Information about target.
			- 'references': Table with information about reference stars close to target.
			- 'avoid': Table with stars close to target which should be avoided in FOV selection.

	.. codeauthor:: Rasmus Handberg <rasmush@phys.au.dk>
	"""

	assert output in ('table', 'json', 'dict'), "Invalid output format"

	# Get API token from config file:
	config = load_config()
	token = config.get('api', 'token', fallback=None)
	if token is None:
		raise Exception("No API token has been defined")

	#
	r = requests.get('https://neo.phys.au.dk/pipeline/reference_stars.php',
		params={'target': target},
		headers={'Authorization': 'Bearer ' + token})
	r.raise_for_status()
	jsn = r.json()

	if output in ('json', 'dict'):
		return jsn

	dict_tables = {}

	tab = Table(
		names=('targetid', 'target_name', 'ra', 'decl', 'catalog_downloaded'),
		dtype=('int32', 'str', 'float64', 'float64', 'bool'),
		rows=[jsn['target']])

	tab['ra'].description = 'Right ascension'
	tab['ra'].unit = u.deg
	tab['decl'].description = 'Declination'
	tab['decl'].unit = u.deg
	dict_tables['target'] = tab


	for table_name in ('references', 'avoid'):
		tab = Table(
		   names=('starid', 'ra', 'decl', 'pm_ra', 'pm_dec', 'gaia_mag', 'gaia_bp_mag', 'gaia_rp_mag', 'gaia_variability', 'H_mag', 'J_mag', 'K_mag', 'g_mag', 'r_mag', 'i_mag', 'z_mag', 'distance'),
		   dtype=('int64', 'float64', 'float64', 'float32', 'float32', 'float32', 'float32', 'float32', 'int32', 'float32', 'float32', 'float32', 'float32', 'float32', 'float32', 'float32', 'float64'),
		   rows=jsn[table_name])

		tab['starid'].description = 'Unique identifier in REFCAT2 catalog'
		tab['ra'].description = 'Right ascension'
		tab['ra'].unit = u.deg
		tab['decl'].description = 'Declination'
		tab['decl'].unit = u.deg
		tab['pm_ra'].description = 'Proper motion in right ascension'
		tab['pm_ra'].unit = u.mas/u.yr
		tab['pm_dec'].description = 'Proper motion in declination'
		tab['pm_dec'].unit = u.mas/u.yr
		tab['distance'].description = 'Distance from object to target'
		tab['distance'].unit = u.deg

		tab['gaia_mag'].description = 'Gaia G magnitude'
		tab['gaia_bp_mag'].description = 'Gaia Bp magnitude'
		tab['gaia_rp_mag'].description = 'Gaia Rp magnitude'
		tab['gaia_variability'].description = 'Gaia variability classification'
		tab['H_mag'].description = '2MASS H magnitude'
		tab['J_mag'].description = '2MASS J magnitude'
		tab['K_mag'].description = '2MASS K magnitude'
		tab['g_mag'].description = 'g magnitude'
		tab['r_mag'].description = 'r magnitude'
		tab['i_mag'].description = 'i magnitude'
		tab['z_mag'].description = 'z magnitude'

		# Add some meta-data to the table as well:
		tab.meta['targetid'] = int(dict_tables['target']['targetid'])

		dict_tables[table_name] = tab

	return dict_tables
