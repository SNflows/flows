#!/usr/bin/env python3
'''
  codeauthor:: Emir Karamehmetoglu <emir.k@phys.au.dk>
'''
import astropy.io.fits as fits
import numpy as np
from astropy.wcs import WCS
import argparse

def main():
    parser = argparse.ArgumentParser(prog='FixTemplateSubHeaderforLCOGT',
                                     description='Attach header and extra dimensions from second \
                                     file to the first, while checking to make sure WCS is identical.\
                                     If not, uses WCS offirst file.')

    parser.add_argument("diff_file", help="Specify the diff file here")
    parser.add_argument("orig_file", help="Specify the original file here")
    parser.add_argument("--overwrite", '-o', type=bool, default=False,
                        help='Overwrite difffile. Default=False outputs new file [diff_file_basename]_hfix.fits')

    args = parser.parse_args()

    hdul_diff = fits.open(args.diff_file)
    wcs_diff = WCS(hdul_diff[0].header)

    hdul_orig = fits.open(args.orig_file)
    wcs_orig = WCS(hdul_orig[0].header)

    if wcs_diff.low_level_wcs == wcs_diff.low_level_wcs:
        new_hdul = fits.HDUList()

        new_hdul.append(hdul_diff[0])
        new_hdul.append(hdul_orig[1])
        new_hdul.append(hdul_orig[2])
        if args.overwrite:
            print('Warning: Overwriting diff file!')
            new_hdul.writeto(args.diff_file,overwrite=True)
        else:
            print('WCS is matched writing to [dill_file]_hfix.fits')
            new_hdul.writeto(args.diff_file.split('.')[:-1][0]+'_hfix.fits',overwrite=True)
    else:
        print('WCS does not match, there may be problems in reducing the file \
               if extension 2 & 3 needs to have correct WCS. \
               Writing to file anyway..')
        new_hdul = fits.HDUList()

        new_hdul.append(hdul_diff[0])
        new_hdul.append(hdul_orig[1])
        new_hdul.append(hdul_orig[2])
        if args.overwrite:
            print('Warning: Overwriting diff file!')
            new_hdul.writeto(args.diff_file,overwrite=True)
        else:
            new_hdul.writeto(args.diff_file.split('.')[:-1][0]+'_hfix.fits',overwrite=True)
main()
