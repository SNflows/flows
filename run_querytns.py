#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Query TNS for new targets and upload to candidate marshal.
https://wis-tns.weizmann.ac.il/
TNS bot apikey must exist in config

.. codeauthor:: Emir Karamehmetoglu <emir.k@phys.au.dk>
.. codeauthor:: Rasmus Handberg <rasmush@phys.au.dk>
"""

import argparse
import logging
from flows import api, load_config
import requests
import json
from collections import OrderedDict
from bs4 import BeautifulSoup
import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.time import Time
from datetime import datetime, timedelta
import time

#--------------------------------------------------------------------
# TNS API FUNCTIONS
# Pre-provided helper functions for the TNS API
# Obtained from https://wis-tns.weizmann.ac.il/content/tns-getting-started
# function for changing data to json format
url_tns_api = "https://wis-tns.weizmann.ac.il/api/get"
def format_to_json(source):
    # change data to json format and return
    parsed = json.loads(source,object_pairs_hook=OrderedDict)
    result = json.dumps(parsed,indent=4)
    return result

# function for search obj
def search(url,json_list):
    try:
        # url for search obj
        search_url = url+'/search'
        # change json_list to json format
        json_file = OrderedDict(json_list)
        # construct the list of (key,value) pairs
        search_data = [('api_key',(None, api_key)),
                     ('data',(None,json.dumps(json_file)))]
        # search obj using request module
        response = requests.post(search_url, files=search_data)
        # return response
        return response
    except Exception as e:
        return [None,'Error message : \n'+str(e)]

# function for get obj
def get_obj_tns(url,json_list):
    try:
        # url for get obj
        get_url = url+'/object'
        # change json_list to json format
        json_file = OrderedDict(json_list)
        # construct the list of (key,value) pairs
        get_data = [('api_key',(None, api_key)),
                     ('data',(None,json.dumps(json_file)))]
        # get obj using request module
        response = requests.post(get_url, files=get_data)
        # return response
        return response
    except Exception as e:
        return [None,'Error message : \n'+str(e)]

# Get SN names from TNS html table
def get_names(soup1):
    names = soup1.find_all('td','cell-name')
    names_list = []
    for name in names:
        if name.text.startswith('SN'):
            names_list.append(name.text.replace(' ',''))
    return names_list

#--------------------------------------------------------------------
if __name__ == '__main__':
    # Parse command line arguments:
    parser = argparse.ArgumentParser(description='Query TNS and upload to candidate marshal')
    parser.add_argument('-d', '--debug', help='Print debug messages.', action='store_true')
    parser.add_argument('-q', '--quiet', help='Only report warnings and errors.', action='store_true')
    parser.add_argument('-z', '--zmax', type=float, help='maximum redshift', default=0.105)
    parser.add_argument('--zmin', type=float, help='minimum redshift', default=0.000000001)
    parser.add_argument('-b', '--days_begin', type=int, help='Discovery day at least X days before today', default=30)
    parser.add_argument('-e', '--days_end', type=int, help='Discovery day at most X days before today', default=3)
    parser.add_argument('-o', '--objtype', type=str, help='TNS objtype int given as comma separed string with no spaces', default='3,104')
    parser.add_argument('-m', '--limit_months', type=int, help='Integer number of months to limit TNS search (for speed). \
                        Should be greater than days_begin', default=2)
    args = parser.parse_args()

    # Set logging level:
    logging_level = logging.INFO
    verbose = True
    if args.quiet:
        logging_level = logging.WARNING
        verbose = False
    elif args.debug:
        logging_level = logging.DEBUG
        verbose = True

    # Setup logging:
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    console = logging.StreamHandler()
    console.setFormatter(formatter)
    logger = logging.getLogger(__name__)
    if not logger.hasHandlers():
        logger.addHandler(console)
    logger.setLevel(logging_level)

    # API key for Bot
    config = load_config()
    api_key = config.get('TNS', 'api_key')
    # PARAMETERS
    days_before_today = args.days_end  # end = now - dbt
    days_to_include = args.days_begin  # begin = now - dti
    z_min = args.zmin
    z_max = args.zmax
    months = args.limit_months  # pre-limit TNS search to candidates reported in the last X months
    if np.floor(days_before_today/30) > months:
        logger.warning('Months limit restricts days_begin, consider increasing limit_months.')
    objtype = args.objtype  # Relevant TNS SN Ia subtypes.

    #CALCULATE CURRENT DATE AND DATE RANGE TO SEARCH
    dt_end = datetime.date(datetime.utcnow()) - timedelta(days=days_before_today)
    dt_begin = datetime.date(datetime.utcnow()) - timedelta(days=days_to_include)
    date_begin = dt_begin.isoformat()
    date_end = dt_end.isoformat()

    #QUERY STRING
    query = "https://wis-tns.weizmann.ac.il/search?&discovered_period_value={}\
&discovered_period_units=months&unclassified_at=0&classified_sne=1\
&include_frb=0&name=&name_like=0&isTNS_AT=all&public=all&ra=&decl=&radius=\
&coords_unit=arcsec&reporting_groupid[]=null&groupid[]=null\
&classifier_groupid[]=null\
&objtype[]={}\
&at_type[]=null\
&date_start[date]={}\
&date_end[date]={}&discovery_mag_min=&discovery_mag_max=\
&internal_name=&discoverer=&classifier=&spectra_count=&redshift_min={}\
&redshift_max={}&hostname=&ext_catid=&ra_range_min=&ra_range_max=\
&decl_range_min=&decl_range_max=&discovery_instrument[]=null\
&classification_instrument[]=null&associated_groups[]=null\
&at_rep_remarks=&class_rep_remarks=&frb_repeat=all&frb_repeater_of_objid=\
&frb_measured_redshift=0&frb_dm_range_min=&frb_dm_range_max=&frb_rm_range_min=\
&frb_rm_range_max=&frb_snr_range_min=&frb_snr_range_max=&frb_flux_range_min=\
&frb_flux_range_max=&num_page=500&display[redshift]=1&display[hostname]=1\
&display[host_redshift]=1&display[source_group_name]=1\
&display[classifying_source_group_name]=1&display[discovering_instrument_name]=0\
&display[classifing_instrument_name]=0&display[programs_name]=0\
&display[internal_name]=1&display[isTNS_AT]=0&display[public]=1\
&display[end_pop_period]=0&display[spectra_count]=1\
&display[discoverymag]=1&display[discmagfilter]=1\
&display[discoverydate]=1&display[discoverer]=1\
&display[remarks]=0&display[sources]=0&display[bibcode]=0\
&display[ext_catalogs]=0".format(months,objtype,date_begin,date_end,z_min,z_max)

    logger.info('Date begin = {}, date_end = {}'.format(date_begin,date_end))
    logger.debug(query)

    # Query TNS
    logger.info('querying TNS for all targets, this may take awhile')
    con = requests.post(query)
    con.raise_for_status()

    logger.info('query successful with status code: {}'.format(con.status_code))
    # Parse output using BS4, find html table
    soup = BeautifulSoup(con.text, 'html.parser')
    tab = soup.find("table",'results-table')
    logger.info('HTML table found')

    if tab is None:
        logger.warn('No HTML table obtained from query!')
        raise AttributeError('No HTML Table Found')

    #Get SN names
    nms = get_names(tab)

    # Remove already existing names using flows api
    included_names = []
    flows_targets = api.get_targets()
    for target in flows_targets:
        included_names.append('SN'+target['target_name'])
    nms = list(set(nms)-set(included_names))
    logger.info('Target names obtained: {}'.format(nms))

    # Query TNS for object info using API, then upload to FLOWS using API.
    sleep_time = 0.25 # seconds to wait in between each query
    for name in nms:
        sn = name.replace('SN','')
        logger.info('querying TNS for: {}'.format(sn))

        # make GET request to TNS via API
        get_obj = [("objname",sn), ("photometry","0"), ("spectra","0")]
        response = get_obj_tns(url_tns_api,get_obj)

        # Parse output
        if None not in response:
            response.raise_for_status()

            parsed = response.json()
            #parsed = json.loads(response.text,object_pairs_hook=OrderedDict)
            logger.info('GET query successful')

            # Extract object info
            reply = parsed['data']['reply']
            name = reply['objname']
            coord = SkyCoord(ra=reply['radeg']*u.deg,dec=reply['decdeg']*u.deg,frame='ICRS')
            redshift = reply['redshift']
            discovery_date = Time(reply['discoverydate'],format='iso',scale='utc')
            discovery_mag = reply['discoverymag']
            host_galaxy = reply['hostname']
            ztf = None
            if 'ZTF' in reply['internal_names']:
                rptnms = reply['internal_names'].replace(' ','').split(',')
                ztf = rptnms['ztf' in rptnms]

            # Try to upload to FLOWS
            try:
                newtargetid = api.add_target(name, coord, redshift=redshift, discovery_date=discovery_date,
                               discovery_mag=discovery_mag, host_galaxy=host_galaxy, ztf=ztf,
                               status='candidate', project='flows')
                logger.info('upload to FLOWS with targetid = {}'.format(newtargetid))
            except ValueError:
                logger.warning('upload to FLOWS failed for: {}'.format(name), exc_info=1)
        else:
            logger.info('GET query not successful for {} \nResponse code: {}'.format(sn,response[1]))
        # Wait before next loop
        time.sleep(sleep_time)
