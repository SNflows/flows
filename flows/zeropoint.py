#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Provides functions for computing the zeropoint and its error 
using bootstrapping with sigma clipping as outlier rejection,
where the sigma is determined by the Chauvenet criteria. Also
allows for arbitrary outlier and fitting functions.

.. codeauthor:: Emir Karamehmetoglu <emir.k@phys.au.dk>
"""

import numpy as np
from astropy.stats import bootstrap
from astropy.modeling import fitting
from scipy.special import erfcinv


#Calculate sigma for sigma clipping using Chauvenet
def sigma_from_Chauvenet(Nsamples):
	'''Calculate sigma according to the Cheuvenet criterion'''
	return erfcinv(1./(2*Nsamples)) * (2.)**(1/2)


def bootstrap_outlier(x,y,yerr, n=500, model='None',fitter='None',
	outlier='None', outlier_kwargs={'sigma':3}, summary='median', error='bootstrap',
	parnames=['intercept'], return_vals=True):
	'''x = catalog mag, y = instrumental mag, yerr = instrumental error
	summary = function for summary statistic, np.nanmedian by default.
	model = Linear1D
	fitter = LinearLSQFitter
	outlier = 'sigma_clip'
	outlier_kwargs, default sigma = 3
	return_vals = False will return dictionary
	Performs bootstrap with replacement and returns model.
	'''
	summary = np.nanmedian if summary == 'median' else summary
	error = np.nanstd if error == 'bootstrap' else error

	#Create index for bootstrapping
	ind = np.arange(len(x))

	#Bootstrap indexes with replacement using astropy
	bootstraps = bootstrap(ind,bootnum=n)
	bootstraps.sort() # sort increasing.
	bootinds = bootstraps.astype(int)

	#Prepare fitter
	fitter_instance = fitting.FittingWithOutlierRemoval(fitter(),outlier, **outlier_kwargs)
	#Fit each bootstrap with model and fitter using outlier rejection at each step.
	#Then obtain summary statistic for each parameter in parnames
	pars = {}
	out = {}
	for parname in parnames:
		pars[parname] = np.ones(len(bootinds), dtype=np.float64)
	for i,bs in enumerate(bootinds):
		#w = np.ones(len(x[bs]), dtype=np.float64) if yerr=='None' else (1.0/yerr[bs])**2
		w = (1.0/yerr[bs])**2
		best_fit, sigma_clipped = fitter_instance(model, x[bs], y[bs], weights=w)
		#obtain parameters of interest
		for parname in parnames:
			pars[parname][i] = best_fit.parameters[np.array(best_fit.param_names) == parname][0]
	if return_vals:
		return [summary(pars[par]) for par in pars]

	for parname in parnames:
		out[parname] = summary(pars[parname])
		out[parname+'_error'] = error(pars[parname])
	return out
