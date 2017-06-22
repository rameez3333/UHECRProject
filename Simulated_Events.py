import numpy as np

import CR_funcs as CR

def gen_iso(N):
	N = int(N)
	dts = ["dec", "theta", "phi"]
	dt = [(d, "f") for d in dts]
	events = np.empty(N, dtype = dt)
	sin_decs = 2 * np.random.rand(N) - 1
	decs = np.arcsin(sin_decs)
	events["dec"] = decs
	events["theta"] = CR.dec2theta(decs)
	events["phi"] = 2 * np.pi * np.random.rand(N)
	return events

def gen_iso_exposure(N, ff):
	N = int(N)
	omega_max = max(CR.omega(-np.pi / 2, ff), CR.omega(np.pi / 2, ff))

	count = 0
	decs = np.empty(N)
	while count < N:
		sin_dec = 2 * np.random.rand() - 1
		dec = np.arcsin(sin_dec)
		if omega_max * np.random.rand() > CR.omega(dec, ff): continue
		decs[count] = dec
		count += 1
	dts = ["dec", "theta", "phi"]
	dt = [(d, "f") for d in dts]
	events = np.empty(N, dtype = dt)

	events["dec"] = decs
	events["theta"] = CR.dec2theta(decs)
	events["phi"] = 2 * np.pi * np.random.rand(N)
	return events

def gen_iso_exposure_one_exp(N, Auger = True):
	dts = ["dec", "theta", "phi"]
	dt = [(d, "f") for d in dts]
	events = np.empty(N, dtype = dt)

	if Auger:	omega_max = CR.omega(-np.pi / 2, True, 1)
	else:		omega_max = CR.omega(np.pi / 2, False, 1)

	count = 0
	decs = np.empty(N)
	while count < N:
		sin_dec = 2 * np.random.rand() - 1
		dec = np.arcsin(sin_dec)
		# ff is irrelevant for one exp
		if omega_max * np.random.rand() > CR.omega(dec, Auger, 1): continue
		decs[count] = dec
		count += 1

	events["dec"] = decs
	events["theta"] = CR.dec2theta(decs)
	events["phi"] = 2 * np.pi * np.random.rand(N)
	return events

