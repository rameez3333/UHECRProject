import numpy as np

import CR_funcs as CR

# energy in EeV, angles in radians
def load_Auger(E_min):
	dts			= ["E",	"RA",	"dec",	"l",	"b"]
	cols		= (3,	4,		5,		6,		7)

	dt = [(d, "f") for d in dts]
	dataf = open("inputs/Auger.txt", "r")
	Auger_events = np.loadtxt(dataf, dtype = dt, usecols = cols)
	dataf.close()
	for angle in dts[1:]:
		Auger_events[angle] = Auger_events[angle] * np.pi / 180

	# apply energy cut
	mask = Auger_events["E"] > E_min
	Auger_events = Auger_events[mask]


	# add a dtype
	Auger_events_ = np.empty(Auger_events.shape, dtype = dt + [("is_Auger", "b")])
	for name in Auger_events.dtype.names:
		Auger_events_[name] = Auger_events[name]
	Auger_events_["is_Auger"] = True

	return Auger_events_

# energy in EeV, angles in radians
def load_TA(E_min, purge_hotspot = False):
	dts			= ["E",	"RA",	"dec",	"l",	"b"]
	cols		= (7,	8,		9,		10,		11)

	dt = [(d, "f") for d in dts]
	dataf = open("inputs/TA.txt", "r")
	TA_events = np.loadtxt(dataf, dtype = dt, usecols = cols)
	dataf.close()
	for angle in dts[1:]:
		TA_events[angle] = TA_events[angle] * np.pi / 180

	# apply energy cut
	mask = TA_events["E"] > E_min
	TA_events = TA_events[mask]

	# put TA hotspot back to bkg
	if purge_hotspot:
		# 19 events within 20 degrees compared to 4.5 estimated
		dec_TA_hotspot, RA_TA_hotspot = (43.2 * np.pi / 180, 146.7 * np.pi / 180)
		count_total = 0
		count_removed = 0
		mask = np.ones(len(TA_events), dtype = bool)
		for i in xrange(len(TA_events)):
			if CR.cos_angle_between(dec_TA_hotspot, RA_TA_hotspot, TA_events[i]["dec"], TA_events[i]["RA"]) > np.cos(20 * np.pi / 180):
				count_total += 1
				if np.random.rand() > 4.5 / 19:
					mask[i] = False
					count_removed += 1
		print "TA hotspot events:", count_total, "  Removed", count_removed, "events"
		TA_events = TA_events[mask]

	# add a dtype
	TA_events_ = np.empty(TA_events.shape, dtype = dt + [("is_Auger", "b")])
	for name in TA_events.dtype.names:
		TA_events_[name] = TA_events[name]
	TA_events_["is_Auger"] = False

	return TA_events_

def load_Auger_TA(E_min_Auger = None, E_min_TA = None, E_scale = None, purge_hotspot = False):
	assert [E_min_Auger, E_min_TA, E_scale].count(None) <= 1

	if E_min_TA == None:	E_min_TA = E_min_Auger * E_scale
	if E_min_Auger == None:	E_min_Auger = E_min_TA / E_scale

	Auger_events = load_Auger(E_min_Auger)
	TA_events = load_TA(E_min_TA, purge_hotspot = purge_hotspot)

	# overlap region, calculate fudge factor
	N_Auger_overlap = sum(1 for Auger_event in Auger_events if Auger_event["dec"] > CR.dec_min)
	N_TA_overlap = sum(1 for TA_event in TA_events if TA_event["dec"] < CR.dec_max)

	return np.append(Auger_events, TA_events), CR.fudge(N_Auger_overlap, N_TA_overlap)

