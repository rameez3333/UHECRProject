import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import chi2
from scipy.integrate import quad
from scipy.special import lpmv
from scipy.optimize import fsolve
from healpy.sphtfunc import Alm

import CR_funcs as CR
#import KMatrix as KM

def realYlm(l, m, theta, phi):
	if abs(m) > l: return 0
	theta = np.array(theta)
	phi = np.array(phi)
	x = np.cos(theta)
	if m == 0:
		return np.sqrt((2 * l + 1) / (4 * np.pi)) * lpmv(0, l, x)
	N = np.sqrt((2 * l + 1) / (4 * np.pi)) * np.sqrt(factorial_ratio(l - abs(m), l + abs(m)))
	P = lpmv(abs(m), l, x)
	if m < 0:
		return np.sqrt(2) * N * P * np.sin(abs(m) * phi)
	return np.sqrt(2) * N * P * np.cos(m * phi)

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
	return np.zeros((int(l_max) + 1) ** 2)

def get_alms_iso(events, l_max):
	alms = np.empty((l_max + 1) ** 2)
	for l in xrange(l_max + 1):
		for m in xrange(-l, l + 1):
			alms[lm2i(l, m)] = np.mean(realYlm(l, m, events["theta"], events["phi"]))
	return alms

def alms_exposure_K(events, l_max, Kinv):
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

	exposure_alms = np.empty((l_max + 1) ** 2)
	for l in xrange(l_max + 1):
		for m in xrange(-l, l + 1):
			alms = realYlm(l, m, CR.dec2theta(events["dec"]), events["phi"]) / omegas
			exposure_alms[lm2i(l, m)] = np.sum(alms)
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
	convert denton's complex alms to healpy's real alms
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
