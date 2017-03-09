import numpy as np
from scipy.stats import chi2

def get_chisq(alms1, alms2, N):
	"""
	alms = [(1,-1),(1,0),(1,1),(2,-2),...]
	N is the number of events that lead to the data sample
	"""
	assert len(alms1) == len(alms2)
	l_max = np.sqrt(len(alms1))
	assert l_max == int(l_max)
	alm_std = 1. / np.sqrt(4 * np.pi * N)
	chisq = 0
	for i in xrange(1, len(alms1)): # no information in ell=0
		chisq += ((alms1[i] - alms2[i]) / alm_std) ** 2
	return chisq

def p_value(chisq, l_max):
	dof = (l_max + 1) ** 2 - 1 # -1 because there is no information in the ell = 0 mode
	return chi2.sf(chisq, dof)

if __name__ == "__main__":
	import Data
	import Ylm
	import Simulated_Events

	# set up properties
	E_min_Auger = 57
	E_scale = 1.13
	purge_hotspot = False
	l_max_max = 10

	# load in data
	Auger_TA_events, ff = Data.load_Auger_TA(E_min_Auger = E_min_Auger, E_scale = E_scale, purge_hotspot = purge_hotspot)
	N = len(Auger_TA_events)

	# simulate events without regard to exposure and with exposure
	Simulated_events = Simulated_Events.gen_iso(len(Auger_TA_events))
	Simulated_events_exposure = Simulated_Events.gen_iso_exposure(len(Auger_TA_events), ff)

	# calculate alms
	Exact_iso_alms = Ylm.get_alms_iso_exact(l_max_max)
	Simulated_alms = Ylm.get_alms_iso(Simulated_events, l_max_max)
	Auger_TA_alms = Ylm.get_alms_exposure(Auger_TA_events, l_max_max, ff)
	Simulated_alms_exposure = Ylm.get_alms_exposure(Simulated_events_exposure, l_max_max, ff)

	# calculate chisq's and pvalues, print in a pretty way
	def out(alms1, alms2, N, l_max_max):
		for l_max in xrange(1, l_max_max + 1):
			i = Ylm.lm2i(l_max, l_max)
			chisq = get_chisq(alms1[:i + 1], alms2[:i + 1], N)
			dof = (l_max + 1) ** 2 - 1 # -1 because there is no information in the ell = 0 mode
			print "l_max = {0:2d}, chisq = {1:6.2f}, dof = {2:3d}, chisq/dof = {3:4.2f}, p = {4:9.7f}".format(l_max, chisq, dof, chisq / dof, p_value(chisq, l_max))

	print "N events =", N

	print "Comparing simulated iso in uniform exposure to iso..."
	out(Simulated_alms, Exact_iso_alms, N, l_max_max)
	print ""

	print "Comparing simulated iso in nonuniform exposure to iso..."
	out(Simulated_alms_exposure, Exact_iso_alms, N, l_max_max)
	print ""

	print "Comparing data to iso..."
	out(Auger_TA_alms, Exact_iso_alms, N, l_max_max)
