import re, logging

import numpy as np
import astropy.units as u

from astropy.table import Table

COMMANDS = {
    'maskstar': ( # mask star for galaxy subtracted psf photometry
        '(\d+(?:\.\d*)?),\s?([+-]?\d+(?:\.\d*)?)', # parse
        lambda ra, dec: (float(ra), float(dec)), # convert
        lambda ra, dec, hdul: 0 <= ra < 360 and -90 <= dec <= 90 # validate
    ),
    'localseq': ( # specifiy hdu name with custom local sequence
        '(.+)', # parse
        lambda lsqhdu: (lsqhdu.upper(),), # convert
        lambda lsqhdu, hdul: lsqhdu in [hdu.name for hdu in hdul], # validate
    ),
    'colorterm': ( # color term for moving references to another system
        '([+-]?\d+(?:\.\d*)?)\s?\((.+)-(.+)\)', # parse
        lambda cterm, A_mag, B_mag: (float(cterm), A_mag, B_mag), # convert
        lambda cterm, A_mag, B_mag, hdul: True, # validate
    )
}

def maskstar(data, wcs, stars, fwhm):

    if not stars:
        return

    data = data.copy()
    if not hasattr(data, 'mask'):
        data = np.ma.array(data, mask=np.zeros_like(data))
    
    X, Y = np.meshgrid(*map(np.arange, data.shape[::-1]))

    for x, y in wcs.all_world2pix(stars, 0):
        i = np.where(((X-x)**2 + (Y-y)**2 < fwhm**2))
        data.mask[i] = True

    return data

def localseq(lsqhdus, hdul):

    if not lsqhdus:
        return

    hdu = hdul[[hdu.name for hdu in hdul].index(lsqhdus[-1][0])]

    references = Table(hdu.data)
    n = len(references)

    references.add_column(np.arange(n) + 1, name='starid', index=0)
    references['pm_ra']  = np.zeros(n) * u.deg / u.yr
    references['pm_dec'] = np.zeros(n) * u.deg / u.yr

    references.meta['localseq'] = hdu.name

    return references

def colorterm(ref_filter, colorterms, references):

    if not colorterms:
        return

    colorterm, A_mag, B_mag = colorterms[-1]

    if not {A_mag, B_mag} <= set(references.colnames):
        missing = ', '.join({A_mag, B_mag} - set(references.colnames))
        logging.warning('%s not in references', missing)
        return

    references[ref_filter] += colorterm * (references[A_mag] - references[B_mag])

    return references

def get_fitscmd(image, command):

    if not command in COMMANDS:
        logging.warning('fitscmd %s unknown', command)
        return

    if not '' in image.head:
        return

    parameters = [
            c[13+len(command):] for c in image.head['']
            if c.startswith('FLOWSCMD: %s' % command)
    ]

    return [
            COMMANDS[command][1](*re.match(COMMANDS[command][0], p).groups())
            for p in parameters
    ]
