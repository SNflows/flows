from flows import api, load_config
from astropy.time import Time
import datetime
from flows import catalogs
from astropy.coordinates import SkyCoord
import astropy.units as u
import argparse
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from astropy.wcs import WCS
import regions
from astroquery.skyview import SkyView
from flows.plots import plot_image
from astropy.visualization import ZScaleInterval
from reproject import reproject_interp

sns.set()

## Manually enter targetid here if running in notebook
## Else code will parse from command line arguments.
targetid = 1349
default_pars = (None, targetid, 0.0, 0.0)  # Rotation, tid, alpha, delta,
skip_prop = True


def is_notebook():
	try:
		__IPYTHON__
		return True
	except NameError:
		return False


def parse():
	parser = argparse.ArgumentParser(description='Calculate Brightest Star')
	parser.add_argument('-t', '--targetid', help="calculate for this targetid",
						type=int, default=None, action='append')
	parser.add_argument('-r', '--rot', help='rotation angle in degrees',
						type=float, default=0.0, action='append')
	parser.add_argument('-a', '--shifta', help='shift alpha in arcsec',
						type=float, default=0.0, action='append')
	parser.add_argument('-d', '--shiftd', help='shift delta in arcsec',
						type=float, default=0.0, action='append')

	return parser.rot, parser.targetid, parser.shifta, parser.shiftd


def get_corners(ra, dec, radius):
	# Get corners for given ra dec and side length
	return [(ra - radius, dec - radius), (ra - radius, dec + radius),
			(ra + radius, dec + radius), (ra + radius, dec - radius)]


def process_corner(corners):
	_points = [u.quantity.Quantity(corner) for corner in corners]
	return np.array(_points)


def default_point_chip1(ra, dec, chip1_center):
	return ra + chip1_center, dec - chip1_center


if __name__ == '__main__':
	rot, tid, shifta, shiftd = default_pars if is_notebook() else parse()
	skip_shift = (np.array((shifta, shiftd)) == 0.0).all()  # skip shift if alpha and delta 0
	# define useful values and get target info and catalog.
	today = Time(datetime.datetime.today())
	config = load_config()
	c = api.get_catalog(tid)
	ref = c['references']
	tar = c['target']
	target_info = api.get_target(tid)
	target_info['skycoord'] = SkyCoord(target_info['ra'] * u.deg, target_info['decl'] * u.deg)
	simbad = catalogs.query_simbad(target_info['skycoord'])
	simbad_coords = simbad[1].apply_space_motion(new_obstime=Time(2000, format='decimalyear'))
	# simbad_coords =  simbad[1] if skip_prop else simbad[1].apply_space_motion(new_obstime=today)
	refcoords = SkyCoord(ref['ra'], ref['decl'])

	ra_tar = tar['ra'].data[0]
	dec_tar = tar['decl'].data[0]

	chip1_center = 3.75 * u.arcmin + 15 * u.arcsecond  # This is the distance to Chip 1 center from Field center.
	field_hw = 7.5 * u.arcminute  # full field height width
	chip1_hw = 3.5 * u.arcminute  # chip1 field height width

	# Get default pointing roughly centered at chip1 center.
	ra, dec = default_point_chip1(ra_tar * u.deg, dec_tar * u.deg, chip1_center)

	# Get corners of Chip 1 and full-field centered, where full field is centered at ra,dec.
	corners = get_corners(ra, dec, field_hw)
	chip1 = get_corners(ra - chip1_center,
						dec + chip1_center, chip1_hw)

	# Process to be easily plottable
	chip1 = process_corner(chip1)
	corners = process_corner(corners)

	out = SkyView.get_images(position='{}, {}'.format(tar['ra'].data[0], tar['decl'].data[0]),
							 survey='2MASS-H', pixels='2500',
							 coordinates='J2000', scaling='Linear', radius=50 * u.arcmin)
	image = out[0][0]

	wcs_H = WCS(image.header)
	wcs_H2 = wcs_H.deepcopy()
	ref_pix = wcs_H.all_world2pix(ref['ra'], ref['decl'], 0)

	# Rotation of image.
	# This is disabled for now, we rotate the region boxes instead!
	if rot:
		rot_angle = rot
	else:
		rot_angle = 0.0
	theta = np.radians(0.0)
	rot_matrix = np.array([[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]])
	wcs_H2.wcs.pc = rot_matrix
	array, footprint = reproject_interp((image.data, wcs_H), wcs_H2, shape_out=image.data.shape)
	wcs_H = wcs_H2

	# Regions
	offset_hawki = 112.5 * u.arcsecond
	if skip_shift:
		offset_alpha, offset_delta = offset_hawki, offset_hawki
		new_center_coord = target_info['skycoord']
	else:
		offset_alpha, offset_delta = offset_hawki + shifta * u.arcsecond, offset_hawki + shiftd * u.arcsecond
		new_center_coord = target_info['skycoord'].spherical_offsets_by(-shifta * u.arcsecond, -shiftd * u.arcsecond)

	chip1_center_true = 3.75 / 2 * u.arcminute + 15 * u.arcsecond  # This is the distance to Chip 1 center from Field center.
	field_hw = 7.5 * u.arcminute  # full field height width
	chip1_hw = 3.5 * u.arcminute  # chip1 field height width

	Hawki_center_coord = target_info['skycoord'].spherical_offsets_by(-offset_alpha, -offset_delta)

	sep = Hawki_center_coord.spherical_offsets_to(target_info['skycoord'])  # unused
	field_angle = rot_angle * u.deg

	chip1_center_coord = Hawki_center_coord.directional_offset_by(
		Hawki_center_coord.position_angle(new_center_coord) + field_angle,
		Hawki_center_coord.separation(new_center_coord) + 15 * u.arcsecond
	)

	chip1_r = regions.RectangleSkyRegion(chip1_center_coord, width=chip1_hw, height=chip1_hw, angle=field_angle)
	Hawki_r = regions.RectangleSkyRegion(Hawki_center_coord, width=field_hw, height=field_hw, angle=field_angle)

	Hawki_pixel_region = Hawki_r.to_pixel(wcs_H)
	chip1_pixel_region = chip1_r.to_pixel(wcs_H)

	mask = ref['H_mag'] <= 14.0
	ref_pix_masked = wcs_H.all_world2pix(ref['ra'][mask], ref['decl'][mask], 0)

	corners_pix = wcs_H.all_world2pix(corners[:, 0] * u.deg, corners[:, 1] * u.deg, 0)
	tar_pix = wcs_H.all_world2pix(ra_tar, dec_tar, 0)
	chip1_pix = wcs_H.all_world2pix(chip1[:, 0], chip1[:, 1], 0)
	simbad_stars_pix = wcs_H.all_world2pix(simbad_coords.ra, simbad_coords.dec, 0)

	# Check bright stars
	simbad_stars = Hawki_r.contains(simbad_coords, wcs_H)
	catalog_stars = Hawki_r.contains(refcoords, wcs_H)

	bright_stars = np.hstack((simbad[0][simbad_stars]['H_mag'].data, ref[catalog_stars]['H_mag'].data))
	if np.ma.isMaskedArray(bright_stars): bright_stars = bright_stars.data[~np.isnan(bright_stars.data)]

	print('Brightest star has H-mag = {0:.1f}'.format(np.round(bright_stars.min(), 1)))

	# Plot
	zscale = ZScaleInterval()
	vmin, vmax = zscale.get_limits(image.data.flat)
	sns.set_theme(style='ticks', font_scale=2.5)
	fig, ax = plt.subplots(figsize=(20, 20), subplot_kw={'projection': wcs_H}, dpi=200)

	Hawki_pixel_region.plot(ax=ax, edgecolor='cyan', linestyle='-.', label='HAWK-I')
	chip1_pixel_region.plot(ax=ax, edgecolor='black', linestyle='-.', label='chip1')
	plot_image(array, ax=ax, cmap='viridis', scale='linear', vmin=vmin, vmax=vmax)

	ref_pix_masked = wcs_H.all_world2pix(ref['ra'][mask], ref['decl'][mask], 0)
	ax.scatter(ref_pix_masked[0], ref_pix_masked[1],
			   facecolors='none', edgecolors='red', zorder=5, alpha=3.0, marker='o', s=200, label='bright stars')

	corners_pix = wcs_H.all_world2pix(corners[:, 0] * u.deg, corners[:, 1] * u.deg, 0)
	tar_pix = wcs_H.all_world2pix(ra_tar, dec_tar, 0)
	ax.scatter(tar_pix[0], tar_pix[1], marker='*', s=450, label='SN', color='orange')

	chip1_pix = wcs_H.all_world2pix(chip1[:, 0], chip1[:, 1], 0)
	simbad_stars = wcs_H.all_world2pix(simbad_coords.ra, simbad_coords.dec, 0)
	ax.scatter(simbad_stars[0], simbad_stars[1],
			   facecolors='none', edgecolors='orange', zorder=5, alpha=3.0, marker='s', s=200, label='simbad stars')

	ax.legend(fontsize=18)
	fig.savefig('{}_HAWKI_H_FC.png'.format(target_info['target_name']))
