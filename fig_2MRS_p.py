import matplotlib.pyplot as plt
import numpy as np

import Data
import Ylm
import Chisq

import warnings

warnings.simplefilter("error")

E_min_Auger = 53
E_scale = 1.13
purge_hotspot = False
l_max = 1
fname = "2MRSEarthmap.txt"

events, ff = Data.load_Auger_TA(E_min_Auger = E_min_Auger, E_scale = E_scale, purge_hotspot = purge_hotspot)
cplx_alms_Auger_TA = Ylm.get_alms_exposure(events, l_max, ff)

real_alms_Auger_TA = Ylm.cplx2real_alms(cplx_alms_Auger_TA)
real_alms_2MRS = Ylm.fname2real_alms(fname, l_max)

#print real_alms_Auger_TA
#print real_alms_2MRS

chisq = Chisq.get_chisq(real_alms_Auger_TA, real_alms_2MRS, len(events))
p = Chisq.p_value(chisq, l_max)
print chisq, p
exit()

chisqs = np.zeros(l_max + 1)
for l in xrange(1, l_max + 1):
	i_max = Ylm.lm2i(l, l)
	chisqs[l] = Chisq.get_chisq(real_alms_Auger_TA[:i_max + 1], real_alms_2MRS[:i_max + 1], len(events))

l_maxs = xrange(1, l_max + 1)

p_value_cumulative = [Chisq.p_value(chisqs[l_max], l_max = l_max) for l_max in l_maxs]
p_value_individual = [Chisq.p_value(chisqs[l_max] - chisqs[l_max - 1], dof = 2 * l_max + 1) for l_max in l_maxs]

for l_max in l_maxs: print l_max, chisqs[l_max] - chisqs[l_max - 1], 2 * l_max + 1, (chisqs[l_max] - chisqs[l_max - 1]) / (2 * l_max + 1), p_value_individual[l_max - 1], chisqs[l_max], (l_max + 1) ** 2 - 1, chisqs[l_max] / ((l_max + 1) ** 2 - 1), p_value_cumulative[l_max - 1]

plt.plot(l_maxs, p_value_cumulative, "-", label = r"${\rm Cumulative}$")
plt.plot(l_maxs, p_value_individual, "--", label = r"${\rm Individual}$")

#print p_value_cumulative
#print p_value_individual

plt.legend(title = r"$\ell:$", loc = 4, fontsize = 14)

plt.xlabel(r"$\ell$")
plt.ylabel(r"$p$")

#plt.yscale("log")

ax = plt.gca().twiny()
ax.set_xlim(plt.gca().get_xlim())
upper_tick_names = np.array([180, 60, 30, 15, 10, 8, 7, 6, 5])
upper_tick_locations = 180. / upper_tick_names
ax.set_xticks(upper_tick_locations)
ax.set_xticklabels([r"$%i$" % upper_tick_name for upper_tick_name in upper_tick_names])
ax.set_xlabel(r"$180/\ell$")

plt.savefig("fig/2MRS_p.eps")

