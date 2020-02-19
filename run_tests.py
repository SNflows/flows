# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 22:42:58 2020

@author: au195407
"""

import matplotlib.pyplot as plt
import numpy as np
from flows.catalogs import download_catalog, get_catalog

#--------------------------------------------------------------------------------------------------
if __name__ == '__main__':

	#download_catalog()

	refs = get_catalog(2)

	target = refs['target']
	ref_ra = np.array([a['ra'] for a in refs['references']])
	ref_dec = np.array([a['decl'] for a in refs['references']])
	ref_jmag = np.array([a['J_mag'] for a in refs['references']])
	avoid_ra = np.array([a['ra'] for a in refs['avoid']])
	avoid_dec = np.array([a['decl'] for a in refs['avoid']])

	plt.figure()
	r = plt.scatter(ref_ra, ref_dec, alpha=0.5, c=ref_jmag)
	plt.scatter(avoid_ra, avoid_dec, marker='x', color='r')
	plt.scatter(target['ra'], target['decl'], marker='+', color='b')
	plt.axis('equal')
	plt.colorbar(r)
	plt.show()