#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
  code author: Emir Karamehmetoglu (emir.k@phys.au.dk)

Query ZTF target information using alerce.
https://alerceapi.readthedocs.io/en/latest/ztf_db.html
'''
from flows import api
import requests
import json
import pandas as pd

#Get target info from flows
cat = api.get_catalog(8, output='json')
tgt = cat['target']

#Check if OID exists in catalog
try:
    oid = tgt['oid']
except KeyError:
    print('ZTF object id not found, querying alerce api with coordinates.')

    #Make json query for Alerce query API
    query = {
        "records_per_pages":20,
        "query_parameters": {
          "coordinates":{
            "ra":  ra,
            "dec": dec,
            "sr": 0.1
        }
        },
        "total":5
    }

    #Run http POST json query to alerce following their API
    url = 'http://ztf.alerce.online/query'
    response = requests.post(url,json=query)
    print("Status code: ", response.status_code)

    #If successful, get objectid of first object
    if response.status_code == 200:
        df = response.json()
        names = list(df['result'].keys())
        oid = df['result'][names[0]]['oid']

#Query ALERCE for detections of object based on oid
url = 'http://ztf.alerce.online/get_detections'
response = requests.post(url,json={'oid':oid})
det = response.json()['result']['detections']
df2 = pd.DataFrame(det)

#Get LC data out and save as CSV files
bands = {}
fid_dict = {1:'g',2:'r',3:'i'}
for fid in df2.fid.unique():
    bands[fid_dict[fid]] = df2[['mjd','magpsf_corr','sigmapsf_corr']][df2.fid == fid]
for band in bands:
    bands[band].reset_index(drop=True,inplace=True)
    bands[band]['jd'] = bands[band]['mjd'] + 2400000.5
    bands[band].columns = ['mjd','mag','err','jd']
    bands[band].to_csv(tgt['target_name'] + '_ztf' + band + '.csv')
