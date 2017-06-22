import numpy as np
import matplotlib.pyplot as plt
import healpy as hp
from astropy import units as u
from astropy.coordinates import SkyCoord

pscoords = np.genfromtxt('PS_1_1_p30_u300_E57.0_1.txt', delimiter='|').transpose()
arr = np.genfromtxt('maps_LH1_1_p30_u300_E57.0_1nside32.txt', delimiter='|')
pickedarr = np.genfromtxt('UHECR_1_1_p30_u300_E57.0_1.txt', delimiter='|').transpose()

pickedcoord = SkyCoord(l = pickedarr[1]*u.radian, b = pickedarr[0]*u.radian, frame='galactic')
pickedra = pickedcoord.icrs.ra.value
pickeddec = pickedcoord.icrs.dec.value   

hp.mollview(arr[2], title='Theta Map')
hp.projscatter(pscoords[0],pscoords[1], lonlat=True, marker='*', color='red')
hp.projscatter(pickedra,pickeddec, lonlat=True, color='black')
plt.show()
