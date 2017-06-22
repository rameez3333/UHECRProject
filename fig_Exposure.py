import matplotlib.pyplot as plt
import numpy as np

import CR_funcs as CR
import Data

E_min_Auger = 57
E_scale = 1.13

tmp, ff = Data.load_Auger_TA(E_min_Auger = E_min_Auger, E_scale = E_scale, purge_hotspot = False)

decs = np.linspace(-np.pi / 2, np.pi / 2, 100)
omegas_auger = CR.omega_one_exp(decs, True, ff)
omegas_TA = CR.omega_one_exp(decs, False, ff)
omegas_both = omegas_auger + omegas_TA

decs *= 180 / np.pi
plt.plot(decs, omegas_both, "k-")
plt.plot(decs, omegas_auger, "k:")
plt.plot(decs, omegas_TA, "k:")

plt.xlabel(r"${\rm Declination}$")
plt.ylabel(r"${\rm Relative\ Exposure}$")

v = list(plt.axis())
v[0] = -90
v[1] = 90
v[2] = 0
plt.axis(v)

plt.xticks([-90, -45, 0, 45, 90])

plt.text(-55, 0.3 * v[3], r"${\rm Auger}$")
plt.text(65, 0.14 * v[3], r"${\rm TA}$")

plt.savefig("fig/Exposure.eps")

