#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Plotting utilities.

.. codeauthor:: Rasmus Handberg <rasmush@phys.au.dk>
"""

import logging
import os
import numpy as np
from bottleneck import allnan, anynan
import matplotlib
from matplotlib.ticker import MaxNLocator
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from astropy.visualization import (PercentileInterval, ImageNormalize, SqrtStretch,
	LogStretch, LinearStretch)

# Change to a non-GUI backend since this
# should be able to run on a cluster:
plt.switch_backend('Agg')

#--------------------------------------------------------------------------------------------------
def plot_image(image, scale='log', origin='lower', xlabel='Pixel Column Number',
	ylabel='Pixel Row Number', make_cbar=False, clabel='Flux ($e^{-}s^{-1}$)', cbar_ticks=None, cbar_ticklabels=None,
	title=None, percentile=95.0, vmin=None, vmax=None, ax=None, cmap=None, offset_axes=None, **kwargs):
	"""
	Utility function to plot a 2D image.

	Parameters:
		image (2d array): Image data.
		scale (str or astropy.visualization.ImageNormalize object, optional): Normalization used to stretch the colormap.
			Options: ``'linear'``, ``'sqrt'``, or ``'log'``. Can also be a `astropy.visualization.ImageNormalize` object.
			Default is ``'log'``.
		origin (str, optional): The origin of the coordinate system.
		xlabel (str, optional): Label for the x-axis.
		ylabel (str, optional): Label for the y-axis.
		make_cbar (boolean, optional): Create colorbar? Default is ``False``.
		clabel (str, optional): Label for the color bar.
		title (str or None, optional): Title for the plot.
		percentile (float, optional): The fraction of pixels to keep in color-trim. The same fraction of pixels is eliminated from both ends. Default=95.
		ax (matplotlib.pyplot.axes, optional): Axes in which to plot. Default (None) is to use current active axes.
		cmap (matplotlib colormap, optional): Colormap to use. Default is the ``Blues`` colormap.
		kwargs (dict, optional): Keyword arguments to be passed to `matplotlib.pyplot.imshow`.
	"""

	if allnan(image):
		logger = logging.getLogger(__name__)
		logger.error("Image is all NaN")
		return None

	# Special treatment for boolean arrays:
	if isinstance(image, np.ndarray) and image.dtype == 'bool':
		if vmin is None: vmin = 0
		if vmax is None: vmax = 1
		if cbar_ticks is None: cbar_ticks = [0, 1]
		if cbar_ticklabels is None: cbar_ticklabels = ['False', 'True']

	#elif not isinstance(image, np.ma.MaskedArray) and anynan(image):
	#	image = np.ma.masked_array(image, mask=np.isnan(image))

	# Calcualte limits of color scaling:
	if vmin is None or vmax is None:
		vmin1, vmax1 = PercentileInterval(percentile).get_limits(image)
		if vmin is None: vmin = vmin1
		if vmax is None: vmax = vmax1

	# Create ImageNormalize object with extracted limits:
	if scale == 'log':
		norm = ImageNormalize(vmin=vmin, vmax=vmax, stretch=LogStretch())
	elif scale == 'linear':
		norm = ImageNormalize(vmin=vmin, vmax=vmax, stretch=LinearStretch())
	elif scale == 'sqrt':
		norm = ImageNormalize(vmin=vmin, vmax=vmax, stretch=SqrtStretch())
	elif isinstance(scale, matplotlib.colors.Normalize) or isinstance(scale, ImageNormalize):
		norm = scale
	else:
		raise ValueError("scale {} is not available.".format(scale))

	if offset_axes:
		extent = (offset_axes[0]-0.5, offset_axes[0] + image.shape[1]-0.5, offset_axes[1]-0.5, offset_axes[1] + image.shape[0]-0.5)
	else:
		extent = (-0.5, image.shape[1]-0.5, -0.5, image.shape[0]-0.5)

	if ax is None:
		ax = plt.gca()

	if cmap is None:
		cmap = plt.cm.Blues
		cmap.set_bad('k')
	elif isinstance(cmap, str):
		cmap = plt.get_cmap(cmap)

	im = ax.imshow(image, origin=origin, norm=norm, extent=extent, cmap=cmap, interpolation='nearest', **kwargs)
	if xlabel is not None: ax.set_xlabel(xlabel)
	if ylabel is not None: ax.set_ylabel(ylabel)
	if title is not None: ax.set_title(title)
	ax.set_xlim([extent[0], extent[1]])
	ax.set_ylim([extent[2], extent[3]])

	if make_cbar:
		fig = ax.figure
		divider = make_axes_locatable(ax)
		if make_cbar == 'top':
			#cax = divider.append_axes('top', size='5%', pad=0.05)
			#cax.xaxis.set_ticks_position("top")
			orientation = 'horizontal'
		elif make_cbar == 'bottom':
			#cax = divider.append_axes('bottom', size='5%', pad=0.2)
			orientation = 'horizontal'
		elif make_cbar == 'left':
			#cax = divider.append_axes('left', size='5%', pad=0.2)
			orientation = 'vertical'
		else:
			#cax = divider.append_axes('right', size='5%', pad=0.05)
			orientation = 'vertical'
		
		cbar = fig.colorbar(im, norm=norm, ax=ax, orientation=orientation)
		cbar.set_label(clabel)
		if cbar_ticks is not None: cbar.set_ticks(cbar_ticks)
		if cbar_ticklabels is not None: cbar.set_ticklabels(cbar_ticklabels)

	# Settings for ticks (to make Mikkel happy):
	ax.xaxis.set_major_locator(MaxNLocator(integer=True))
	ax.xaxis.set_minor_locator(MaxNLocator(integer=True))
	ax.yaxis.set_major_locator(MaxNLocator(integer=True))
	ax.yaxis.set_minor_locator(MaxNLocator(integer=True))
	ax.tick_params(direction='out', which='both', pad=5)
	ax.xaxis.tick_bottom()
	#ax.set_aspect(aspect)

	return im

