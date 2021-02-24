from .reference_stars import reference_stars
reference_stars = type('reference_stars', (dict,), {
    'target' : {
        'targetid' : -1,
        'target_name' : None,
        'target_status' : None,
        'ra' : None,
        'decl' : None,
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
    },
    'references' : {
        'starid' : -1,
        'ra' : None,
        'decl' : None,
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
    }
})(reference_stars)
reference_stars.avoid = reference_stars.references

from .datafiles import datafiles
datafiles = type('datafiles', (dict,), {
    'image' : {
        'fileid' : None,
        'path' : None,
        'targetid' : None,
        'site' : None,
        'filesize' : None,
        'filehash' : None,
        'inserted' : '1970-01-01 00:00:00.0',
        'lastmodified' : '1970-01-01 00:00:00.0',
        'photfilter' : None,
        'obstime' : None,
        'exptime' : None,
        'version' : None,
        'archive_path' : None,
        'target_name' : None,
        'template' : None,
        'diffimg' : None
    },
    'diffimg' : {
        'fileid' : None,
        'path' : None,
        'filehash' : None,
        'filesize' : None
    }
})(datafiles)

from .sites import sites
sites = type('sites', (dict,), {
    'site' : {
        'siteid' : -1,
        'sitename' : None,
        'longitude': None,
        'latitude' : None,
        'elevation' : None,
        'site_keyword' : None,
    }
})(sites)
