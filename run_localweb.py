import os, io

import numpy as np
import matplotlib.pyplot as plt

from flask import Flask, session; session = dict() # XXX
from base64 import b64encode as b64
from astropy.io import ascii

from flows import api, load_config

app = Flask(__name__)
#app.secret_key = os.urandom(24)

@app.route('/')
def index():
    s = str()
    target_by_name = {target['target_name']: target for target in api.get_targets()}
    local_targets = os.listdir(load_config()['photometry']['archive_local'])
    for target_name in sorted(target_by_name):
        targetid = target_by_name[target_name]['targetid']
        s += f'<a href=\'/{targetid}\'>{target_name}</a> \n' if target_name in local_targets else f'{target_name} \n'
    return s

def lightcurve(targetid):
    data = dict()
    target = {target['targetid']: target for target in api.get_targets()}[targetid]
    output = load_config()['photometry']['output']
    if not os.path.exists(output + '/%s' % (target['target_name'])): return ''
    for fileid in os.listdir(output + '/%s' % (target['target_name'])):
        f = output + '/%s/%s/photometry.ecsv' % (target['target_name'], fileid)
        if not os.path.exists(f): continue
        table = ascii.read(f)
        filt, mjd = table.meta['photfilter'], table.meta['obstime-bmjd']
        if not filt in data: data[filt] = []
        mag, err, starid = table[0]['mag'], table[0]['mag_error'], table[0]['starid']
        data[filt].append((fileid, mjd, mag, err, starid))
    plt.figure(figsize=(20, 10))
    for filt in data:
        fileid, mjd, mag, err, starid = map(np.array, zip(*data[filt]))
        _ = plt.errorbar([], [], yerr=[], ls='', marker='o', label=filt)
        plt.errorbar(mjd[starid==-1], mag[starid==-1], yerr=err[starid==-1], ls='', marker='o', color=_.lines[0].get_color())
        plt.errorbar(mjd[starid==0], mag[starid==0], yerr=err[starid==0], ls='', marker='s', color=_.lines[0].get_color())
        for i in range(len(fileid)): plt.text(mjd[i], mag[i], str(int(fileid[i])), fontsize=8)
    plt.title(target['target_name'])
    plt.xlabel('MJD'); plt.ylabel('Magnitude')
    plt.legend()
    plt.tight_layout()
    ymin, ymax = min([min(list(zip(*d))[2]) for d in data.values()]), max([max(list(zip(*d))[2]) for d in data.values()])
    plt.ylim(ymin - (ymax-ymin)*0.05, ymax + (ymax-ymin)*0.05)
    plt.gca().invert_yaxis()
    png = io.BytesIO()
    plt.savefig(png, format='png', dpi=100)
    plt.close()
    return b64(png.getvalue()).decode('utf-8')

@app.route('/<int:targetid>')
def target(targetid):
    s = '<img width=\'100%\' src=\'data:image/png;base64, {}\' />'.format(lightcurve(targetid))
    target = {target['targetid']: target for target in api.get_targets()}[targetid]
    archive = os.listdir(load_config()['photometry']['archive_local'] + '/%s' % target['target_name'])
    output = load_config()['photometry']['output']
    fileids = list(map(str, api.get_datafiles(targetid, 'all'))) # XXX
    for fileid in set(fileids) - set(session.keys()): session[fileid] = api.get_datafile(int(fileid))
    s += '<table border=1>'
    for fileid in fileids:
        s += '<tr><td>{fileid}</td><td>{path}</td><td>{obstime}</td><td>{photfilter}</td>'.format(**session[fileid])
        if session[fileid]['path'].rsplit('/',1)[1] in archive:
            if not session[fileid]['path'].rsplit('/',1)[1][:-8] + '.png' in archive: s += '<td>IMG</td>'
            else: s += f'<td><a href=\'/{targetid}/{fileid}/img\'>IMG</a></td>'
        else: s += '<td></td>'
        if not os.path.exists(output + '/%s/%0.5d/photometry.log' % (target['target_name'], int(fileid))): s += '<td>LOG</td>'
        else: s += f'<td><a href=\'/{targetid}/{fileid}/log\'>LOG</a></td>'
        if not os.path.exists(output + '/%s/%0.5d/photometry.ecsv' % (target['target_name'], int(fileid))): s += '<td>PHOT</td>'
        else: s += f'<td><a href=\'/{targetid}/{fileid}/phot\'>PHOT</a></td>'
        s += '</tr>'
    s += '</table>'
    return s

@app.route('/<int:targetid>/<int:fileid>/img')
def image(targetid, fileid):
    s = str()
    target = {target['targetid']: target for target in api.get_targets()}[targetid]
    archive = load_config()['photometry']['archive_local']
    output = load_config()['photometry']['output']
    img = '%s/%s.png' % (archive, session[str(fileid)]['path'][:-8])
    with open(img, 'rb') as fd: img = b64(fd.read()).decode('utf-8')
    s += f'<img width=\'25%\' src=\'data:image/png;base64, {img}\' />'
    if not os.path.exists(output + '/%s/%0.5d' % (target['target_name'], fileid)): return s
    for f in sorted(os.listdir(output + '/%s/%0.5d' % (target['target_name'], fileid)))[::-1]:
        if not '.' in f or not f.rsplit('.', 1)[1] == 'png': continue
        img = output + '/%s/%0.5d/%s' % (target['target_name'], fileid, f)
        with open(img, 'rb') as fd: img = b64(fd.read()).decode('utf-8')
        s += f'<img width=\'25%\' src=\'data:image/png;base64, {img}\' />'
    return s

@app.route('/<int:targetid>/<int:fileid>/phot')
def photometry(targetid, fileid):
    target = {target['targetid']: target for target in api.get_targets()}[targetid]
    output = load_config()['photometry']['output']
    phot = output + '/%s/%0.5d/photometry.ecsv' % (target['target_name'], fileid)
    table = ascii.read(phot)
    s = str(table.meta) + '<br /><br /><table>'
    s += '<tr><td><b>' + '</b></td><td><b>'.join(table.colnames) + '</b></td></tr>'
    for row in table: s += '<tr><td>' + '</td><td>'.join(map(str, row)) + '</td></tr>'
    s += '</table>'
    return s


@app.route('/<int:targetid>/<int:fileid>/log')
def log(targetid, fileid):
    target = {target['targetid']: target for target in api.get_targets()}[targetid]
    output = load_config()['photometry']['output']
    log = output + '/%s/%0.5d/photometry.log' % (target['target_name'], fileid)
    with open(log, 'r') as fd: return fd.read().replace('\n', '<br />')

if __name__ == '__main__':
    app.run(debug=True)
