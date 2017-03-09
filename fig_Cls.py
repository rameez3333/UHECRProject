import matplotlib.pyplot as plt
import numpy as np

import Data
import Ylm

E_min_Auger = 57
E_scale = 1.13
E_min_TA = E_min_Auger * E_scale

l_max = 25

Auger_TA_events, ff = Data.load_Auger_TA(E_min_Auger = E_min_Auger, E_scale = E_scale, purge_hotspot = False)
Auger_TA_events_purged, ff_purged = Data.load_Auger_TA(E_min_Auger = E_min_Auger, E_scale = E_scale, purge_hotspot = True)

alms = Ylm.get_alms_exposure(Auger_TA_events, l_max, ff)
alms_purged = Ylm.get_alms_exposure(Auger_TA_events_purged, l_max, ff_purged)
Cls = Ylm.get_Cls(alms)
Cls_purged = Ylm.get_Cls(alms_purged)

plt.plot(xrange(1, l_max + 1), Cls[1:], "b-") # a00 is always 1/sqrt(4pi)
#plt.plot(xrange(1, l_max + 1), Cls_purged[1:], "r-", label = r"${\rm TA\ Hotspot\ Subtracted}$")

plt.xlabel(r"$\ell$")
plt.ylabel(r"$C_\ell$")

v = list(plt.axis())
v = [0, l_max + 1, 0, v[3]]
plt.axis(v)

ells = np.linspace(1, l_max + 1, 100)
p = 0.95
Cl_95_ranges = np.array([Ylm.Cl_region(l, len(Auger_TA_events), p) for l in ells])
plt.fill_between(ells, Cl_95_ranges[:, 0], Cl_95_ranges[:, 1], lw = 0, facecolor = (0.8, 0.8, 0.8), label = r"$%i\%%{\rm\ CL\ isotropic}$" % (100 * p))

plt.legend(loc = 4, numpoints = 1, frameon = False, fontsize = 14)

plt.savefig("fig/Cls.eps")

