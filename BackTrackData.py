# backtracks the data through the GMF

import numpy as np
import matplotlib.pyplot as plt

import crpropa as crp
import healpy as hp
from astropy.coordinates import SkyCoord
from astropy import units as u

import Data
import CR_funcs as CR

def draw_map_lines():
	hp.visufunc.graticule(dpar = 30, dmer = 30)
	hp.visufunc.projtext(np.pi / 2, +0.99 * np.pi, r"$+180^\circ\ $", ha = "right")
	hp.visufunc.projtext(np.pi / 2, -0.99 * np.pi, r"$\ -180^\circ$", ha = "left")
	plt.text(0.98, 0.02, r"${\rm Galactic}$", transform = plt.gca().transAxes, ha = "right", va = "bottom", fontsize = 12)

# parameters
E_min_Auger = 53
E_scale = 1.13
E_min_TA = E_min_Auger * E_scale
purge_hotspot = False
nside = 128
Auger_angular_sigma = np.deg2rad(0.9) # 1 sigma angular uncertainty in degrees for Auger
TA_angular_sigma = np.deg2rad(1.35) # 1 sigma (I guess) angular uncertainty in degrees for TA, the average of 1 and 1.7
gal_size = 20 # size of the sphere outside the galaxy to propagate to in kpc
seed = 1888
rng = crp.Random(seed)
nrepeat_backtrack = 1e3 # number of back propagated cr's from each data point
# either calculate or read in from file:
read_in_galaxy_map = True
read_in_deflection_map = True

# load data
events, ff = Data.load_Auger_TA(E_min_Auger = E_min_Auger, E_scale = E_scale, purge_hotspot = purge_hotspot)
Auger_events = Data.load_Auger(E_min_Auger)
TA_events = Data.load_TA(E_min_TA, purge_hotspot = purge_hotspot)

# draw map at earth in equatorial coordinates
# initialize maps
Auger_map = np.zeros(hp.nside2npix(nside))
TA_map = np.zeros(hp.nside2npix(nside))
earth_map = np.zeros(hp.nside2npix(nside))

# fill in maps
for event in Auger_events:	Auger_map[hp.pixelfunc.ang2pix(nside, np.pi / 2 - event["b"], event["l"])] += 1. / CR.omega_one_exp(event["dec"], True, ff)
for event in TA_events:		TA_map[hp.pixelfunc.ang2pix(nside, np.pi / 2 - event["b"], event["l"])] += 1. / CR.omega_one_exp(event["dec"], False, ff)

# smooth maps
earth_map += hp.sphtfunc.smoothing(Auger_map, sigma = Auger_angular_sigma, verbose = False)
earth_map += hp.sphtfunc.smoothing(TA_map, sigma = TA_angular_sigma, verbose = False)

# maximum direction
theta_gal, phi_gal = hp.pix2ang(nside, np.argmax(earth_map))
coord = SkyCoord(b = (np.pi / 2 - theta_gal) * u.rad, l = phi_gal * u.rad, frame = "galactic")
print "Maximum is at", coord.icrs

# plot
hp.visufunc.mollview(earth_map, cbar = False, title = r"${\rm Data\ at\ Earth}$")

# draw eq plane
n = 100
decs = np.zeros(n)
RAs = np.linspace(0, 2 * np.pi, n)
coords = SkyCoord(dec = decs * u.rad, ra = RAs * u.rad, frame = "icrs")
bs = coords.galactic.b.rad
ls = coords.galactic.l.rad
mp = hp.projector.MollweideProj()
xs, ys = mp.ang2xy(np.pi / 2 - bs, ls)

# remove the break
for i in xrange(1, n):
	if abs(xs[i] - xs[i - 1]) > 1: break
xs = np.roll(xs, -i)
ys = np.roll(ys, -i)
#plt.plot(xs, ys, "k-")

draw_map_lines()
plt.savefig("fig/earth_data.eps", bbox_inches = "tight")
plt.clf()

# crpropa backtracking
# set up GMF
B = crp.JF12Field()
B.randomStriated(seed)
B.randomTurbulent(seed)

# set up simulation
sim = crp.ModuleList()
sim.add(crp.PropagationCK(B, 1e-3, 0.1 * crp.parsec, 100 * crp.parsec))
obs = crp.Observer()
obs.add(crp.ObserverLargeSphere(crp.Vector3d(0), gal_size * crp.kpc))
sim.add(obs)

position = crp.Vector3d(-8.5, 0, 0) * crp.kpc # start at the earth

pid = crp.nucleusId(1, 1) # proton
pid *= -1 # antiparticle, for backtracking

if read_in_galaxy_map:
	galaxy_map = np.genfromtxt("maps/galaxy_map.txt")
else:
	galaxy_map = np.zeros(hp.nside2npix(nside))
	for i in xrange(len(events)):
		if i % 10 == 0: print i, len(events)
		event = events[i]
		if event["is_Auger"]:
			mean_energy = event["E"] * crp.EeV
			omega = CR.omega_one_exp(event["dec"], True, ff)
			sigma = Auger_angular_sigma
		else:
			mean_energy = event["E"] * crp.EeV / E_scale
			omega = CR.omega_one_exp(event["dec"], False, ff)
			sigma = TA_angular_sigma
		sigma_energy = 0.2 * mean_energy # see fig 19 in 1407.3214 for Auger, bottom of page 5 of 1404.5890 for TA
		mean_direction = crp.Vector3d()
		mean_direction.setRThetaPhi(1, np.pi / 2 - event["b"], 1.0 * event["l"])

		for i in xrange(int(nrepeat_backtrack)):
			energy = rng.randNorm(mean_energy, sigma_energy)
			initial_direction = rng.randVectorAroundMean(mean_direction, sigma)
			c = crp.Candidate(crp.ParticleState(pid, energy, position, initial_direction))
			sim.run(c)
			final_direction = c.current.getDirection()
			galaxy_map[hp.pixelfunc.ang2pix(nside, final_direction.getTheta(), final_direction.getPhi())] += 1. / (nrepeat_backtrack * omega)

	np.savetxt("maps/galaxy_map.txt", galaxy_map)
	
hp.visufunc.mollview(galaxy_map, cbar = False, title = r"${\rm Data\ Outside\ the\ Galaxy}$")
draw_map_lines()
plt.savefig("fig/galaxy_data.eps", bbox_inches = "tight")
plt.clf()

# generate a map of mean deflections
nside_deflection = 128
if read_in_deflection_map:
	deflection_map = np.genfromtxt("maps/deflection_map.txt")
else:
	deflection_map = np.zeros(hp.nside2npix(nside_deflection)) # map of the magnitude of the deflections positioned at the earth, to verify that the GC is in the GC

	initial_direction = crp.Vector3d()
	for i in xrange(hp.nside2npix(nside_deflection)):
		if i % 1e3 == 0: print i, hp.nside2npix(nside_deflection)

		theta, phi = hp.pixelfunc.pix2ang(nside_deflection, i)
		initial_direction.setRThetaPhi(1, theta, phi)

		mean_deflection_angle = 0
		for j in xrange(len(events)):
			event = events[j]
			if event["is_Auger"]:
				energy = event["E"] * crp.EeV
			else:
				energy = event["E"] * crp.EeV / E_scale
			c = crp.Candidate(crp.ParticleState(pid, energy, position, initial_direction))
			sim.run(c)
			final_direction = c.current.getDirection()
			mean_deflection_angle += np.rad2deg(initial_direction.getAngleTo(final_direction))
		mean_deflection_angle /= len(events)
		deflection_map[i] = mean_deflection_angle

	np.savetxt("maps/deflection_map.txt", deflection_map)

#deflection_map = hp.pixelfunc.ud_grade(deflection_map, nside_out = nside) # upgrade for smoothing
#deflection_map = hp.sphtfunc.smoothing(deflection_map, sigma = np.deg2rad(2), verbose = False) # add in smoothing
hp.visufunc.mollview(deflection_map, cbar = True, title = r"${\rm Mean\ Deflection\ Angle\ (Degrees)}$")
draw_map_lines()
plt.savefig("fig/deflection_angle.eps", bbox_inches = "tight")
plt.clf()

