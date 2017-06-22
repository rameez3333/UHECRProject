import healpy as hp
import numpy as np
import Data
import CR_funcs as CR

# load in data and calculate fudge factor
E_min_Auger = 57
E_scale = 1.13
events, ff = Data.load_Auger_TA(E_min_Auger = E_min_Auger, E_scale = E_scale)

def logLikelihood(m):
	"""
	takes a healpy map m and returns the likelihood that the UHECR data agrees with it
	assumes that the map is normalized to:
		int dOmega m(Omega) = 1
	"""
	n_side = hp.pixelfunc.npix2nside(len(m))

	logL = 0
	for i in xrange(len(events)):
		dec = events["dec"][i]
		theta = np.pi / 2 - dec
		RA = events["RA"][i]
		phi = RA
		pixel = hp.pixelfunc.ang2pix(n_side, theta, phi)
		logL += np.log(m[pixel] / CR.omega(dec, ff)[0])
	return logL

def test():
	alms = np.array([1+0j, 0.1+0j, 0j])
	n_side = 32
	m = hp.sphtfunc.alm2map(alms, n_side)
	print logLikelihood(m)

if __name__ == "__main__": test()

