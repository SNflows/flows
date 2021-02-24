import numpy as np

from flask import Flask, request
from web import sites, datafiles, reference_stars

import json

app = Flask(__name__)

@app.route('/api/sites.php')
def api_sites():
    if 'siteid' in request.args:
        siteid = int(request.args['siteid'])
        return sites[siteid]
    for site in sites:
        sites[site] = {**sites.site, **sites[site]}
    return json.dumps(list(sites.values()))

@app.route('/api/datafiles.php')
def api_datafiles():
    if 'fileid' in request.args:
        fileid = int(request.args['fileid'])
        for f in (f for target in datafiles for f in datafiles[target]):
            if f['fileid'] != fileid:
                continue
            f.update({**datafiles.image, **f})
            if not f['diffimg'] is None:
                f['diffimg'] = {**datafiles.diffimg, **f['diffimg']}
            return f
    elif 'targetid' in request.args:
        targetid = int(request.args['targetid'])
        return str([f['fileid'] for f in datafiles[targetid]])

@app.route('/api/reference_stars.php')
def api_reference_stars():
    targetid = int(request.args['target'])
    reference_stars[targetid]['target'] = {
        **reference_stars.target, **reference_stars[targetid]['target']
    }
    for reference in reference_stars[targetid]['references']:
        reference.update({**reference_stars.references, **reference})
    for avoid in reference_stars[targetid]['avoid']:
        avoid.update({**reference_stars.avoid, **avoid})
    return reference_stars[targetid]

if __name__ == '__main__':
    app.run(debug=True)
