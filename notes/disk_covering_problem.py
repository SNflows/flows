#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""


.. codeauthor:: Rasmus Handberg <rasmush@phys.au.dk>
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord, SkyOffsetFrame
from matplotlib.patches import Circle
import astropy.units as u

#--------------------------------------------------------------------------------------------------
if __name__ == '__main__':

	coo_centre = SkyCoord(ra=0, dec=0, unit='deg', frame='icrs')

	radius = 24.0/60.0

	#aframe = SkyOffsetFrame(origin=coo_centre)
	#c = coo_centre.transform_to(aframe)
	#print(c)

	fig, ax = plt.subplots()
	ax.plot(coo_centre.ra.deg, coo_centre.dec.deg, 'rx')

	ax.add_artist(Circle([coo_centre.ra.deg, coo_centre.dec.deg], radius=radius, ec='r', fc=None, fill=False))
	ax.add_artist(Circle([coo_centre.ra.deg, coo_centre.dec.deg], radius=0.5*radius, ec='b', fc=None, fill=False))

	for n in range(6):
		new = SkyCoord(
			ra=coo_centre.ra.deg + 0.8 * radius * np.cos(n*60*np.pi/180),
			dec=coo_centre.dec.deg + 0.8 * radius * np.sin(n*60*np.pi/180),
			unit='deg', frame='icrs')

		ax.plot(new.ra.deg, new.dec.deg, 'bx')
		ax.add_artist(Circle([new.ra.deg, new.dec.deg], radius=0.5*radius, ec='b', fc=None, fill=False))


	plt.axis('equal')
	#ax.set_xlim(coo_centre.ra.deg + radius * np.array([-2, 2]))
	#ax.set_ylim(coo_centre.dec.deg +radius * np.array([-2, 2]))
	plt.show()
