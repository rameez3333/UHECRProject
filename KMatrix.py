import numpy as np
from scipy.integrate import quad
from scipy.special import lpmv

import CR_funcs as CR
import Ylm

def K_element(l1, m1, l2, m2, ff, inverse = False):
	# inverse calculates int Ylm Ylm/omega while no inverse calculates in Ylm Ylm * omega
	if m1 != m2: return 0

	# with phi part 'integrated out'
	def realYlm(l, m, theta):
		if abs(m) > l: return 0
		x = np.cos(theta)
		if m == 0:
			return np.sqrt((2 * l + 1.) / 2) * lpmv(0, l, x) # time sqrt(2pi) because phi integral is 2pi
		N = np.sqrt(2 * l + 1) / 2 * np.sqrt(CR.factorial_ratio(l - abs(m), l + abs(m))) # times sqrt(pi) because phi integral is pi
		P = lpmv(abs(m), l, x)
		return np.sqrt(2) * N * P

	if inverse:
		return quad(lambda theta: realYlm(l1, m1, theta) * realYlm(l2, m2, theta) * np.sin(theta) / omega(CR.dec2theta(theta), ff), 0, np.pi, epsrel = 1e-4, limit = 100)[0]
	else:
		return quad(lambda theta: realYlm(l1, m1, theta) * realYlm(l2, m2, theta) * np.sin(theta) * omega(CR.dec2theta(theta), ff), 0, np.pi, epsrel = 1e-4, limit = 100)[0]

def K_matrix(l_max, fudge_factor, inverse = False):
	K = np.matrix(np.zeros(((l_max + 1) ** 2, (l_max + 1) ** 2)))
	for i in xrange((l_max + 1) ** 2):
		l1, m1 = Ylm.i2lm(i)
		for j in xrange(i, (l_max + 1) ** 2):
			l2, m2 = Ylm.i2lm(j)
			element = K_element(l1, m1, l2, m2, fudge_factor, inverse)
			K[i, j] = element
			if j > i: K[j, i] = element # symmetric matrix
	return K

def almKs(blms, Kinv):
	# takes the 1/omega inverse
	l_max_blms = np.sqrt(len(blms)) - 1
	assert l_max_blms == int(l_max_blms)
	assert Kinv.shape[0] == Kinv.shape[1] # is square
	l_max_K = np.sqrt(Kinv.shape[0]) - 1
	assert l_max_K == int(l_max_K)
	assert l_max_K >= l_max_blms
	Kinv = Kinv[:len(blms), :len(blms)] # prune down to size
	almKs = Ylm.get_alms_iso_exact(l_max_blms) # initialize to zero

	for i in xrange(len(almKs)):
		for j in xrange(len(almKs)):
			almKs[i] += Kinv[i, j] * blms[j] # K is symmetric, so Kinv is symmetric
	return almKs / (almKs[0] * np.sqrt(4 * np.pi)) # correct the normalization: a00 = 1/sqrt(4pi)

