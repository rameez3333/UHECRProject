import matplotlib.pyplot as plt
import numpy as np

import Data
import Ylm

E_min_Auger = 57
E_scale = 1.13
E_min_TA = E_min_Auger * E_scale

l_max = 25

Auger_TA_events, ff = Data.load_Auger_TA(E_min_Auger = E_min_Auger, E_scale = E_scale, purge_hotspot = False)

alms = Ylm.get_alms_exposure(Auger_TA_events, l_max, ff)

xs = np.empty(len(alms))
for i in xrange(len(alms)):
	l, m = Ylm.i2lm(i)
	if l == 0:
		xs[i] = 0
	else:
		xs[i] = l + 0.2 * m / l

plt.plot(xs[1:], alms[1:], "b+") # a00 is always 1/sqrt(4pi)

plt.xlabel(r"$\ell+0.2\frac m\ell$")
plt.ylabel(r"$a_{\ell m}$")

v = list(plt.axis())
y_max = max(abs(v[2]), abs(v[3]))
v = [0, l_max + 1, -y_max, y_max]
plt.axis(v)

alm_std = 1. / np.sqrt(4 * np.pi * len(Auger_TA_events))
# 95%Cl => 1.96 sigma
scale = 1.96
plt.fill_between(v[:2], [-scale * alm_std, -scale * alm_std], [scale * alm_std, scale * alm_std], lw = 0, facecolor = (0.8, 0.8, 0.8), label = r"$95\%{\rm\ CL\ isotropic}$")

plt.plot(v[:2], [0, 0], "k-")

plt.legend(loc = 1, frameon = False, fontsize = 14)

plt.savefig("fig/alms.eps")

