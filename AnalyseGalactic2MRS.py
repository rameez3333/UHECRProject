import numpy as np
import matplotlib.pyplot as plt
import healpy as hp
from glob import glob
from astropy import units as u
from astropy.coordinates import SkyCoord

twomrsprojectedmap = 0.046*np.genfromtxt('2MRSProjectedmapSi__0.3N_0.3C_0.25Li_0.1He_0.05H_Be-9.txt') + 0.281*np.genfromtxt('2MRSProjectedmapN__0.3N_0.3C_0.25Li_0.1He_0.05H_Be-9.txt') + 0.673*np.genfromtxt('2MRSProjectedmapHe__0.4He_0.6H_Be-9.txt')

#twomrsprojectedmap = np.genfromtxt('2MRSProjectedMapPureProton.txt')

twomrsprojectedmap = hp.pixelfunc.ud_grade(twomrsprojectedmap, 128, pess=True, power=1.)
twomrsprojectedmap = twomrsprojectedmap/twomrsprojectedmap.sum()
#smoothing by 2 degrees
twomrsprojectedmap = hp.smoothing(twomrsprojectedmap, sigma = np.deg2rad(3.))
twomrsprojectedmap = twomrsprojectedmap-1.*np.min(twomrsprojectedmap)
twomrsprojectedmap = twomrsprojectedmap/twomrsprojectedmap.sum()

#hp.mollview(twomrsprojectedmap, title='2MRS projected map, 2 degree smoothing')
#plt.show()
#plt.savefig('2MRSprojectedmap.png')

def scattomap(dec,ra, nside=32):
    hmap = np.bincount(hp.ang2pix(nside, np.deg2rad(90.-dec), np.deg2rad(ra)), minlength=hp.nside2npix(nside))
    return hmap    

def weightedscattomap(dec,ra, weights, nside=32):
    hmap = np.bincount(hp.ang2pix(nside, np.deg2rad(90.-dec), np.deg2rad(ra)), weights=weights, minlength=hp.nside2npix(nside))
    return hmap

fwcs = ['Outputs/56_26/MinEn53.0_N10_SpecI4.3_Seq60?GalBack_Manual.txt','Outputs/14_7/MinEn53.0_N10_SpecI4.3_Seq60?GalBack_Manual.txt', 'Outputs/12_6/MinEn53.0_N10_SpecI4.3_Seq70?GalBack_Manual.txt', 'Outputs/7_3/MinEn53.0_N10_SpecI4.3_Seq70?GalBack_Manual.txt', 'Outputs/4_2/MinEn53.0_N10_SpecI4.3_Seq70?GalBack_Manual.txt', 'Outputs/1_1/MinEn53.0_N10_SpecI4.3_Seq70?GalBack_Manual.txt']
cweights = [0., 0.3, 0.3, 0.25, 0.1, 0.05]

#cweights = [1.0]
#fwcs = ['Outputs/1_1/MinEn53.0_N10_SpecI4.3_Seq70?GalBack_Manual.txt']

earthmap = np.zeros(hp.nside2npix(128))

for fwc, cweight in zip(fwcs, cweights):
    if not cweight:
        continue    
    flist = sorted(glob(fwc))

    print len(flist), 'files found'

    arr = np.genfromtxt(flist[0], delimiter="|").transpose()
    print len(arr[0]), 'events found'


    for fname in flist[1:]:
        print fname
        now = np.genfromtxt(fname, delimiter="|").transpose()
        print len(now[0]), 'events found'
        arr = np.append(arr, now , axis=1)

    coord = SkyCoord(l = arr[1]*u.radian, b = arr[0]*u.radian, frame='galactic')
    ra = coord.icrs.ra.value
    dec = coord.icrs.dec.value

    defcoord = SkyCoord(l = arr[4]*u.radian, b = arr[3]*u.radian, frame='galactic')
    defra = defcoord.icrs.ra.value
    defdec = defcoord.icrs.dec.value

    defmap = scattomap(defdec, defra, nside=128)
    #hp.mollview(defmap, title='Backpropagated Map at Galactic border, isotropic at Earth')
    #plt.show()
    #plt.savefig('BackProp.png')

    defpixs = hp.ang2pix(128, np.deg2rad(90.-defdec), np.deg2rad(defra))
    weights = twomrsprojectedmap[defpixs]
    contribmap = weightedscattomap(dec, ra, weights, 128)
    contribmap = contribmap/np.sum(contribmap)
    earthmap = earthmap + cweight*contribmap



#np.savetxt('2MRSEarthmapFe_0.01Fe_0.3N_0.3C__0.25Li_0.09He_0.05H.txt', earthmap, delimiter="|")

np.savetxt('2MRSEarthmapAugerBestFit__0.3N_0.3C__0.25Li_0.1He_0.05H.txt', earthmap, delimiter="|")

hp.mollview(earthmap, title='Map at Earth, backpropagated weighted by 2MRS sources')
plt.show()
plt.savefig('2MRSEarthmapAugerBestFit__0.3N_0.3C__0.25Li_0.1He_0.05H.png')
