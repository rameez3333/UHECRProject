import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad

a_0_Auger = -35.2 * np.pi / 180
theta_m_Auger = 80 * np.pi / 180
a_0_TA = 39.3 * np.pi / 180
theta_m_TA = 55 * np.pi / 180
dec_min = a_0_TA - theta_m_TA
dec_max = a_0_Auger + theta_m_Auger

def omega_one_exp(dec, Auger, fudge_factor):
	if Auger:
		a_0 = a_0_Auger
		theta_m = theta_m_Auger
	else:
		a_0 = a_0_TA
		theta_m = theta_m_TA

	dec = np.array(dec)
	if dec.shape == ():
		dec = np.array([dec])

	xi = (np.cos(theta_m) - np.sin(a_0) * np.sin(dec)) / (np.cos(a_0) * np.cos(dec))

	alpha_m = np.empty(len(dec))
	for i in xrange(len(xi)):
		if xi[i] < -1:
			alpha_m[i] = np.pi
		elif xi[i] > 1:
			alpha_m[i] = 0
		else:
			alpha_m[i] = np.arccos(xi[i])

	if Auger:
		correction = fudge_factor
	else:
		correction = 1
	return correction * (np.cos(a_0) * np.cos(dec) * np.sin(alpha_m) + alpha_m * np.sin(a_0) * np.sin(dec))

def omega(dec, ff):
	return omega_one_exp(dec, True, ff) + omega_one_exp(dec, False, ff)
	
def fudge(N_Auger, N_TA):
	# use the uncorrected exposures here
	int_TA = quad(lambda dec: np.cos(dec) * omega_one_exp(dec, False, 1), dec_min, dec_max)[0]
	int_Auger = quad(lambda dec: np.cos(dec) * omega_one_exp(dec, True, 1), dec_min, dec_max)[0]
	return N_Auger * int_TA / (N_TA * int_Auger)

def cos_angle_between(dec1, RA1, dec2, RA2):
	theta1 = dec2theta(dec1)
	phi1 = RA1
	theta2 = dec2theta(dec2)
	phi2 = RA2
	return np.cos(theta1) * np.cos(theta2) + np.cos(phi1 - phi2) * np.sin(theta1) * np.sin(theta2)

# same as vice versa
def dec2theta(dec):
	return np.pi / 2 - dec

def r_angles(N, ff = 1, detector = "None"):
	ret = np.empty(N, dtype = [("dec", "f"), ("RA", "f")])
	ret["RA"] = 2 * np.pi * np.random.rand(N)

	detector = detector.lower()
	Augers = ["auger", "pao", "a", "p"]
	TAs = ["ta", "t"]
	Boths = ["both", "b"]

	if detector in Augers:
		omega_max = omega_one_exp(-np.pi / 2, True, ff)
	if detector in TAs:
		omega_max = omega_one_exp(np.pi / 2, False, ff)
	if detector in Boths:
		omega_max = omega(-np.pi / 2, True, ff)

	if detector in [None, "none", "n"]:
		decs = np.arcsin(2 * np.random.rand(N) - 1)
	else:
		count = 0
		while count < N:
			dec = np.arcsin(2 * np.random.rand() - 1)
			if detector in Augers:
				if omega_max * np.random.rand() > omega_one_exp(dec, True, ff):		continue
			if detector in TAs:
				if omega_max * np.random.rand() > omega_one_exp(dec, False, ff):	continue
			if detector in Boths:
				if omega_max * np.random.rand() > omega(dec, ff):					continue
			ret["dec"][count] = dec
			count += 1
	return ret

