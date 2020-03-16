# -*- coding: utf-8 -*-
"""

.. codeauthor:: Rasmus Handberg <rasmush@phys.au.dk>
"""

import os.path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.dates import DateFormatter
import astropy.units as u
from astropy.time import Time
from datetime import datetime
from astropy.coordinates import SkyCoord, AltAz, get_sun, get_moon
from astropy.visualization import quantity_support
from . import api

#--------------------------------------------------------------------------------------------------
def visibility(target, siteid, date=None, output=None):
	"""
	Create visibility plot.

	Parameters:
		target (str or int):
		siteid (int): Identifier of site.
		date (str, optional): Date for which to create visibility plot.
			Default it to use the current date.
		output (str, optional): Path to file or directory where to place the plot.
			If not given, the plot will be created in memory, and can be shown on screen.

	.. codeauthor:: Rasmus Handberg <rasmush@phys.au.dk>
	"""

	if date is None:
		date = datetime.utcnow().strftime('%Y-%m-%d')

	cat = api.get_catalog(target, output='json')
	tgt = cat['target']
	targetid = tgt['targetid']

	site = api.get_site(siteid)

	# Coordinates of object:
	obj = SkyCoord(ra=tgt['ra'], dec=tgt['decl'], unit='deg', frame='icrs')

	# Observatory:
	observatory = site['EarthLocation']
	utcoffset = (site['longitude']*u.deg/(360*u.deg)) * 24*u.hour

	# Create timestamps to calculate for:
	midnight = Time(date + ' 00:00:00', scale='utc') - utcoffset
	delta_midnight = np.linspace(-12, 12, 1000)*u.hour
	times = midnight + delta_midnight

	# AltAz frame:
	AltAzFrame = AltAz(obstime=times, location=observatory)

	# Object:
	altaz_obj = obj.transform_to(AltAzFrame)

	# The Sun and Moon:
	altaz_sun = get_sun(times).transform_to(AltAzFrame)
	altaz_moon = get_moon(times).transform_to(AltAzFrame)

	sundown_astro = (altaz_sun.alt < -6*u.deg)
	if np.any(sundown_astro):
		min_time = np.min(times[sundown_astro]) - 2*u.hour
		max_time = np.max(times[sundown_astro]) + 2*u.hour
	else:
		min_time = times[0]
		max_time = times[-1]

	quantity_support()
	fig, ax = plt.subplots(1, 1, figsize=(15,9), squeeze=True)
	plt.grid(ls=':', lw=0.5)
	ax.plot(times.datetime, altaz_sun.alt, color='y', label='Sun')
	ax.plot(times.datetime, altaz_moon.alt, color=[0.75]*3, ls='--', label='Moon')
	objsc = ax.scatter(times.datetime, altaz_obj.alt, c=altaz_obj.az, label=tgt['target_name'], lw=0, s=8, cmap='twilight')
	ax.fill_between(times.datetime, 0*u.deg, 90*u.deg, altaz_sun.alt < -0*u.deg, color='0.5', zorder=0) # , label='Night'
	ax.fill_between(times.datetime, 0*u.deg, 90*u.deg, altaz_sun.alt < -18*u.deg, color='k', zorder=0) # , label='Astronomical Night'

	plt.colorbar(objsc, ax=ax, pad=0.01).set_label('Azimuth [deg]')
	ax.legend(loc='upper left')
	ax.minorticks_on()
	ax.set_xlim(min_time.datetime, max_time.datetime)
	ax.set_ylim(0*u.deg, 90*u.deg)
	ax.set_title("%s - %s - %s" % (str(tgt['target_name']), date, site['sitename']), fontsize=14)
	plt.xlabel('Time [UTC]', fontsize=14)
	plt.ylabel('Altitude [deg]', fontsize=16)
	fig.autofmt_xdate()

	formatter = DateFormatter('%d/%m %H:%M')
	ax.xaxis.set_major_formatter(formatter)

	if output:
		datestr = datetime.strptime(date, 'Y-m-d').strftime('Ymd')
		if os.path.isdir(output):
			plotpath = os.path.join(output, "visibility_%s_target%04d_site%02d.png" % (targetid, datestr, siteid))
		else:
			plotpath = output
		fig.savefig(plotpath, bbox_inches='tight', transparent=True)
		plt.close(fig)
		return plotpath

	plt.show()
	return ax
