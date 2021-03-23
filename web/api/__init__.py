import os, json

DIR = os.path.dirname(os.path.abspath(__file__))

with open("%s/targets.json" % DIR, 'r') as fd:
    targets = json.load(fd)

for target in targets:
    targets[target] = {**{
        'target_status' : None,
        'redshift' : None,
        'redshift_error' : None,
        'discovery_mag' : None,
        'catalog_downloaded' : None,
        'pointing_model_created' : '1970-01-01 00:00:00.0',
        'inserted' : '1970-01-01 00:00:00.0',
        'discovery_date' : '1970-01-01 00:00:00.0',
        'project' : None,
        'host_galaxy' : None,
        'ztf_id' : None
    }, **targets[target]}

with open("%s/catalogs.json" % DIR, 'r') as fd:
    catalogs = json.load(fd)

for target in catalogs:
    catalogs[target] = {
        'references': [
            {**{
                'pm_ra' : 0,
                'pm_dec' : 0,
                'gaia_mag' : None,
                'gaia_bp_mag' : None,
                'gaia_rp_mag' : None,
                'J_mag' : None,
                'H_mag' : None,
                'K_mag' : None,
                'g_mag' : None,
                'r_mag' : None,
                'i_mag' : None,
                'z_mag' : None,
                'gaia_variability' : 0,
                'V_mag' : None,
                'B_mag' : None,
                'u_mag' : None,
                'distance' : None
            }, **reference}
        for reference in catalogs[target]]
    }

with open("%s/datafiles.json" % DIR, 'r') as fd:
    datafiles = json.load(fd)

for target in datafiles:
    for i, datafile in enumerate(datafiles[target]):
        if "diffimg" in datafile and datafile["diffimg"]:
            datafile["diffimg"] = {**{
                'filehash' : None,
                'filesize' : None
            }, **datafile["diffimg"]}
        datafiles[target][i] = {**{
            'site' : None,
            'filesize' : None,
            'filehash' : None,
            'inserted' : '1970-01-01 00:00:00.0',
            'lastmodified' : '1970-01-01 00:00:00.0',
            'obstime' : None,
            'exptime' : None,
            'version' : None,
            'archive_path' : None,
            'template' : None,
            'diffimg' : None
        }, **datafile}

with open("%s/sites.json" % DIR, 'r') as fd:
    sites = json.load(fd)

for i, site in enumerate(sites):
    sites[i] = {**{
        'siteid' : -1,
        'sitename' : None,
        'longitude': None,
        'latitude' : None,
        'elevation' : None,
        'site_keyword' : None,
    }, **site}
