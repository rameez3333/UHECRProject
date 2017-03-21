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

