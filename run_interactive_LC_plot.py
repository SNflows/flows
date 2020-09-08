#! /usr/bin/env python3
import os
import glob
import numpy as np
import warnings
#Filter matplotlib depreciation warnings
warnings.filterwarnings("ignore", module='matplotlib')
import matplotlib
import matplotlib.pyplot as plt
from astropy.table import Table
from flows import api, load_config
import argparse
import mplcursors
import seaborn as sns
matplotlib.use('Qt5Agg')


if __name__ == '__main__':
    #Parser:
    parser = argparse.ArgumentParser(description='Plot photometry for target')
    parser.add_argument('--targetid','-t', type=int, required=True, help='Target id: -t <ID>')
    parser.add_argument('--fileid','-i', nargs='*', type=int, default=None,
                        help='[Optional] Specific file ids within target separed by spaces: -i <ID> <ID> <ID>')
    parser.add_argument('--filters','-f',type=str, default=None,
                        help='[Optional] List of comma delimited filters: -f u,b,g,r \n if not provided will use all')
    parser.add_argument('--offset','-jd', type=float, default=2458800.0)
    args = parser.parse_args()

    config = load_config()

    #To use when only plotting some filters
    if args.filters is not None:
        usefilts = args.filters.split(',')

    #To use when only plotting some fileids
    #Parse input fileids:
    if args.fileid is not None:
        # Plot the specified fileid:
        fileids = args.fileid
    else:
        fileids = []
    if len(fileids) > 1:
        #not implemented yet
        pass

    TargetID = args.targetid
    datafiles = api.get_datafiles(TargetID,filt='all')

    #Change to directory, raise if it does not exist
    workdir_root = config.get('photometry', 'output', fallback='.')
    snname = api.get_datafile(datafiles[0])['target_name']
    sndir = os.path.join(workdir_root,snname)
    try:
        os.chdir(sndir)
    except:
        print('No such directory as',sndir)
        raise

    #Get list of photometry files
    phot_files = glob.glob('*/*.ecsv')

    #Load all data into astropytable
    Phot = Table(names=['jd','mag','mag_err','filter','tel','sub','fileid'],dtype=[None,None,None,'S64','S64',bool,'S64'])

    for file in phot_files:
        fileid = file.split('/')[0]
        df = api.get_datafile(fileid) # get fileid info
        jd, filt, tel = df['obstime'] + 2400000.5,df['photfilter'],df['site']

        #get phot of diff image
        AT = Table.read(file)
        AT.add_index('starid')

        try:
            if -1 in AT['starid']:
                mag,mag_err = AT.loc[-1]['mag'],AT.loc[-1]['mag_error']
                sub = True
            elif 0 in AT['starid']:
                print('No subtraction found for:',file,'in filter',filt)
                mag,mag_err = AT.loc[0]['mag'],AT.loc[0]['mag_error']
                sub = False
            else:
                print('No object phot found, skipping: \n',file)
                continue
            Phot.add_row((jd,mag,mag_err,filt,tel,sub,fileid))
        except:
            print('file:',file)
            raise

    #TODO: Use filters argument here
    filters = ['gp','rp','B','V','ip']
    Phots = dict(zip(filters,[Phot[Phot['filter'] == f] for f in filters]))

    #Plot
    sns.set(style='ticks')
    fig,ax = plt.subplots(figsize=(6.4,4),dpi=130)
    fig.subplots_adjust(top=0.95,left=0.1,bottom=0.1,right=0.97)

    shifts = dict(zip(filters,0.0 * np.arange(len(filters))))
    cps = sns.color_palette()
    colors = dict(zip(filters,(cps[2],cps[3],cps[0],cps[-1],cps[1])))

    offset = args.offset
    for filt in filters:
        ax.errorbar(Phots[filt]['jd'] - offset,Phots[filt]['mag'] + shifts[filt],Phots[filt]['mag_err'],
                    marker='s',linestyle='None',label=filt,color=colors[filt])

    ax.invert_yaxis()
    ax.legend()
    ax.set_xlabel('JD-' + str(offset),fontsize=16)
    ax.set_ylabel('App. Mag',fontsize=16)
    ax.set_title(snname)

    mplcursors.cursor(ax).connect(
        "add", lambda sel: sel.annotation.set_text(Phots[str(sel.artist.get_label())]['fileid'][sel.target.index]))
    plt.show(block=True)
