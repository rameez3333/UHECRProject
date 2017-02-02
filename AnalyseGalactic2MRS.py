import numpy as np
import matplotlib.pyplot as plt
import healpy as hp
from glob import glob
from astropy import units as u
from astropy.coordinates import SkyCoord

twomrsprojectedmap = np.genfromtxt('2MRSProjectedmap.txt')

twomrsprojectedmap = hp.pixelfunc.ud_grade(twomrsprojectedmap, 128, pess=True, power=1.)
twomrsprojectedmap = twomrsprojectedmap/twomrsprojectedmap.sum()
#smoothing by 2 degrees
twomrsprojectedmap = hp.smoothing(twomrsprojectedmap, sigma = np.deg2rad(2.))
twomrsprojectedmap = twomrsprojectedmap-1.*np.min(twomrsprojectedmap)
twomrsprojectedmap = twomrsprojectedmap/twomrsprojectedmap.sum()

hp.mollview(twomrsprojectedmap, title='2MRS projected map, 2 degree smoothing')
plt.show()
plt.savefig('2MRSprojectedmap.png')

flist = sorted(glob('Outputs/1_1/MinEn53.0_N10_SpecI4.3_Seq60?GalBack_Manual.txt'))

print len(flist), 'files found'

arr = np.genfromtxt(flist[0], delimiter="|").transpose()
print len(arr[0]), 'events found'


for fname in flist[1:]:
    print fname
    now = np.genfromtxt(fname, delimiter="|").transpose()
    print len(now[0]), 'events found'
    arr = np.append(arr, now , axis=1)
    
    
def scattomap(dec,ra, nside=32):
    hmap = np.bincount(hp.ang2pix(nside, np.deg2rad(90.-dec), np.deg2rad(ra)), minlength=hp.nside2npix(nside))
    return hmap    

def weightedscattomap(dec,ra, weights, nside=32):
    hmap = np.bincount(hp.ang2pix(nside, np.deg2rad(90.-dec), np.deg2rad(ra)), weights=weights, minlength=hp.nside2npix(nside))
    return hmap    


coord = SkyCoord(l = arr[1]*u.radian, b = arr[0]*u.radian, frame='galactic')
ra = coord.icrs.ra.value
dec = coord.icrs.dec.value

defcoord = SkyCoord(l = arr[4]*u.radian, b = arr[3]*u.radian, frame='galactic')
defra = defcoord.icrs.ra.value
defdec = defcoord.icrs.dec.value

defmap = scattomap(defdec, defra, nside=128)
hp.mollview(defmap, title='Backpropagated Map at Galactic border, isotropic at Earth')
plt.show()
plt.savefig('BackProp.png')

defpixs = hp.ang2pix(128, np.deg2rad(90.-defdec), np.deg2rad(defra))

weights = twomrsprojectedmap[defpixs]

earthmap = weightedscattomap(dec, ra, weights, 128)

hp.mollview(earthmap, title='Map at Earth, backpropagated weighted by 2MRS sources')
plt.show()
plt.savefig('Earthweighted.png')
