import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import chi2
from scipy.integrate import quad
from scipy.special import lpmv, sph_harm
from scipy.optimize import fsolve
from healpy.sphtfunc import Alm
import healpy as hp

import CR_funcs as CR
#import KMatrix as KM

def cplxYlm(l, m, theta, phi):
	# phi -> -phi relates scipy with healpy
	return sph_harm(m, l, -phi, theta)

def factorial_ratio(top, bottom):
	if top == bottom:
		return 1
	if bottom > top:
		return 1. / factorial_ratio(bottom, top)
	return top * factorial_ratio(top - 1, bottom)

def lm2i(l, m):
	return l ** 2 + m + l

def i2lm(i):
	l = int(np.sqrt(i))
	m = i - l - l ** 2
	return l, m

def get_alms_iso_exact(l_max):
	assert l_max == int(l_max)
	return np.zeros(Alm.getsize(l_max), dtype = complex)

def get_alms_iso(events, l_max):
	cplx_alms = np.zeros(Alm.getsize(l_max), dtype = complex)
	for l in xrange(l_max + 1):
		for m in xrange(l + 1):
			cplx_alms[Alm.getidx(l_max, l, m)] = np.mean(cplxYlm(l, m, events["theta"], events["phi"]))
	return cplx_alms

def get_alms_exposure_K(events, l_max, Kinv):
	return -1
	# requiring Kinv is slow
	# either K_matrix(inverse = True) or np.linalg.inv(K) gives similar results
	blms = np.empty((l_max + 1) ** 2) # the uncorrected alms
	for l in xrange(l_max + 1):
		for m in xrange(-l, l + 1):
			blms[lm2i(l, m)] = np.sum(realYlm(l, m, CR.dec2theta(events["dec"]), events["phi"])) / len(events)
	alms = KM.almKs(blms, Kinv)
	return alms

def get_alms_exposure(events, l_max, ff):
	# this is the prefered way
	# 1702.07209 eq 8
	# same as Sommers
	omegas = CR.omega(events["dec"], ff)

	# check field names in events
	phi_field = "phi"
	if phi_field not in events.dtype.names:
		names = list(events.dtype.names)
		idx = names.index("RA")
		names[idx] = phi_field
		events.dtype.names = names

	exposure_alms = np.zeros(Alm.getsize(l_max), dtype = complex)
	for l in xrange(l_max + 1):
		for m in xrange(l + 1):
			cplx_alms = cplxYlm(l, m, CR.dec2theta(events["dec"]), events["phi"]) / omegas
			exposure_alms[Alm.getidx(l_max, l, m)] = np.sum(cplx_alms)
	return exposure_alms / (exposure_alms[0] * np.sqrt(4 * np.pi)) # normalize the alms so that a00 = 1/sqrt(4pi)

def get_Cls(alms):
	l_max = np.sqrt(len(alms)) - 1
	assert l_max == int(l_max)
	l_max = int(l_max)
	Cls = np.zeros(l_max + 1)

	for i in xrange(len(alms)):
		l = int(np.sqrt(i))
		Cls[l] += alms[i] ** 2 / (2 * l + 1)
	return Cls

def Cl_std(l, N):
	return (1. / (4 * np.pi * N)) * np.sqrt(((N - 1.) / N) * (2. / (2. * l + 1.)))

def Cl_region(l, N, p):
	"""
	N is the number of events
	p is the probability region, 0.68, 0.95, etc.
	"""
	dof = 2 * l + 1
	def func(x):
		out = [chi2.cdf(x[1], dof) - chi2.cdf(x[0], dof) - p]
		out.append(chi2.pdf(x[0], dof) - chi2.pdf(x[1], dof))
		return out
	init = [0.5 * dof, 1.5 * dof]
	if l < 2:
		init[0] = 0.1
	if l > 30:
		init[0] = 0.6 * dof
	return fsolve(func, init) / (dof * 4 * np.pi * N)

def cplx2real_alms(cplx_alms):
	"""
	convert healpy's complex alms to denton's real alms
	"""
	l_max = Alm.getlmax(len(cplx_alms))
	real_alms = np.zeros((l_max + 1) ** 2)

	for i in xrange(len(cplx_alms)):
		l, m = Alm.getlm(l_max, i)
		if m == 0:
			real_alms[lm2i(l, m)] = cplx_alms[i].real
		else:
			real_alms[lm2i(l, m)] = np.sqrt(2) * cplx_alms[i].real
			real_alms[lm2i(l, -m)] = -(-1) ** m * np.sqrt(2) * cplx_alms[i].imag
	return real_alms

def real2cplx_alms(real_alms):
	"""
	convert denton's real alms to healpy's complex alms
	"""
	l_max = np.sqrt(len(real_alms)) - 1
	assert l_max == int(l_max)
	l_max = int(l_max)

	cplx_alms = np.zeros(Alm.getsize(l_max), dtype = complex)
	for i in xrange(len(real_alms)):
		l, m = i2lm(i)
		if m == 0:
			cplx_alms[Alm.getidx(l_max, l, m)] = real_alms[i]
		if m < 0:
			cplx_alms[Alm.getidx(l_max, l, abs(m))] += -(-1) ** m * 1j * real_alms[i] / np.sqrt(2)
		if m > 0:
			cplx_alms[Alm.getidx(l_max, l, m)] += real_alms[i] / np.sqrt(2)
	return cplx_alms

def fname2real_alms(fname, l_max):
	the_map = np.genfromtxt(fname)
	the_map = hp.pixelfunc.ud_grade(the_map, 128, pess = True, power = 1.)
	the_map /= the_map.sum()
	the_map = hp.smoothing(the_map, sigma = np.deg2rad(2), verbose = False)
	the_map -= 1. * np.min(the_map)
	the_map /= the_map.sum()

#	the_map[hp.pixelfunc.ang2pix(128, np.pi/2 - 0.222, 3.26)] += 0.1 # virgo cluster in equatorial coordinates
#	the_map[hp.pixelfunc.ang2pix(128, 1.23, 5.4)] += 0.1 # cen a in galactic coordinates

	cplx_alms = hp.sphtfunc.map2alm(the_map, l_max)

	cplx_alms /= cplx_alms[Alm.getidx(l_max, 0, 0)]
	cplx_alms /= np.sqrt(4 * np.pi)
	
#	hp.visufunc.mollview(hp.sphtfunc.alm2map(cplx_alms, 128))
#	plt.show()
	
	real_alms = cplx2real_alms(cplx_alms)
	return real_alms

