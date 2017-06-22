import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import itertools
from scipy.interpolate import InterpolatedUnivariateSpline
from astropy.cosmology import Planck15 as cosmo
import math


ethresh = 53.0

def AZtoNCode(A, Z):
    return 1e9 + Z*1e4 + A*1e1

def CodeToAZ(code):
    return int(str(code)[4:6]), int(str(code)[7:9])

inv_comoving_distance = InterpolatedUnivariateSpline(cosmo.comoving_distance(np.linspace(0, 0.1, 100)).value, np.linspace(0, 0.1, 100))

def scattomap(dec,ra, nside=32):
    hmap = np.bincount(hp.ang2pix(nside, np.deg2rad(90.-dec), np.deg2rad(ra)), minlength=hp.nside2npix(nside))
    return hmap

twomrsarr = np.genfromtxt('2MRS/catalog/2mrs_1175_done.dat', skip_header=10, usecols=(1,2,24)).transpose()
twomrstype = np.genfromtxt('2MRS/catalog/2mrs_1175_done.dat', skip_header=10, usecols=(22), dtype="|S5").transpose()
twomrsmag = np.genfromtxt('2MRS/catalog/2mrs_1175_done.dat', skip_header=10, usecols=(8,9,10)).transpose()
print len(twomrsarr[0]), 'objects loaded from the 2MRS catalog'
z = twomrsarr[2]/299792.458


propagated = np.genfromtxt('Outputs/14_7/MinEn53.0MaxEn500.0_InSpecIndex0.9_StartDist300.0_N100000_Seq0.txt')

startindices = np.where(propagated.transpose()[0]==0.)[0]
weights = np.insert(startindices, len(startindices), len(propagated))
weights  = weights[1:] - startindices
counts = weights
weights = weights*1.0
weights = 1./weights
orig = np.asarray(map(CodeToAZ, propagated.transpose()[3]))
det = np.asarray(map(CodeToAZ, propagated.transpose()[1]))

origz, origa = orig.transpose()[1], orig.transpose()[0]
detz, deta = det.transpose()[1], det.transpose()[0]

propagated = propagated.transpose()

propagated = np.vstack([propagated[0], deta, detz, propagated[2], propagated[4]]).transpose()
weights = np.repeat(weights, counts)


weights = weights[propagated.transpose()[3]>ethresh]
propagated=propagated[propagated.transpose()[3]>ethresh]

#hist = np.histogramdd(sample = propagated, weights=weights)



#lengths=[]

#ethresh=100.

#counter=0
#for i in xrange(0, len(startindices)-1):
    #print i
    #try:
        #if np.where(propagated[startindices[i]:startindices[i+1]].transpose()[1]==AZtoNCode(A2, Z2))[0][-1] >= np.where(propagated[startindices[i]:startindices[i+1]].transpose()[2]>ethresh)[0][-1]:
            #print i, np.where(propagated[startindices[i]:startindices[i+1]].transpose()[1]==AZtoNCode(A2, Z2))[0][-1], np.where(propagated[startindices[i]:startindices[i+1]].transpose()[2]>ethresh)[0][-1]
            #lengths.append(propagated[startindices[i]:startindices[i+1]].transpose()[0][np.where(propagated[startindices[i]:startindices[i+1]].transpose()[2]>ethresh)[0][-1]])
            #counter = counter +1
    #except:
        #print 'Well, fuck'

#print counter        



#verify sample is E^-2       
#plt.hist(np.log10(propagated.transpose()[4]), weights=weights*propagated.transpose()[4]*propagated.transpose()[4])
        

#generate a,z histogram

histaz = np.histogram2d(propagated.transpose()[1], propagated.transpose()[2], weights=weights, bins=[27,56], range=[[0,27],[1,57]])
zbins = histaz[1]
abins = histaz[2]
histazobtained = histaz[0]/np.sum(histaz[0])

zindices = np.digitize(propagated.transpose()[1], zbins, right=True)
aindices = np.digitize(propagated.transpose()[2], abins, right=True)

#repeat the process for other starting nuclei also 



histazdesired=histazobtained*0. #just for the shape

#Manually set individual fractions
histazdesired[np.digitize([1], zbins, right=True)[0]][np.digitize([1], abins, right=True)[0]] = 0.05 #40% proton
histazdesired[np.digitize([2], zbins, right=True)[0]][np.digitize([4], abins, right=True)[0]] = 0.09 #30% helium
histazdesired[np.digitize([3], zbins, right=True)[0]][np.digitize([7], abins, right=True)[0]] = 0.25 #30% lithium
histazdesired[np.digitize([6], zbins, right=True)[0]][np.digitize([12], abins, right=True)[0]] = 0.3 #30% Carbon
histazdesired[np.digitize([7], zbins, right=True)[0]][np.digitize([14], abins, right=True)[0]] = 0.3 #30% Nitropgen

#histazdesired[np.digitize([26], zbins, right=True)[0]][np.digitize([56], abins, right=True)[0]] = 0.01# 1% iron

histscaler = histazdesired/histazobtained

compositionweights = histscaler[zindices,aindices]

binwidth = 5.

distbins = np.arange(0, 250., binwidth)

totweights = weights*compositionweights

detenhistob = np.histogram(np.log10(propagated.transpose()[3]), bins =40, weights = totweights, range=[np.min(np.log10(propagated.transpose()[3])), np.max(np.log10(propagated.transpose()[3]))+0.03])
detenhistbins = detenhistob[1]
detenhistob=detenhistob[0]                                                                                                              
detenhistindices = np.digitize(np.log10(propagated.transpose()[3]), bins = detenhistbins, right=True)

inversedetenhistob = 1./detenhistob

#flatten the distribution
detenweights = inversedetenhistob[detenhistindices]
totweights = totweights*detenweights                                                                                             

#now take the spectral index back down to -4.3
totweights = totweights*np.power(propagated.transpose()[3], -4.3)
totweights[np.isinf(totweights)] = 0.
totweights[np.isnan(totweights)] = 0.
#as the plot in the following line shows, the spectrum at this stage still seems to be quite hard. Lets 
histdist = np.histogram(propagated.transpose()[0], weights=totweights, bins = distbins)
plt.hist(propagated.transpose()[0], weights=totweights, bins = distbins)
plt.xlabel('Distance (MPc)')
plt.ylabel('Events')
plt.yscale('log')
plt.savefig('DistancedistroN__0.3N_0.3C_0.25Li_0.1He_0.05H_Be-9.png')
plt.show()
#abins = range(0, 56)
#zbins = range(0, 26)

bins = histdist[1]
histdist=histdist[0]
avmap = np.zeros(hp.nside2npix(32))
bincount=0
print histdist
histdist[0]=histdist[1]
for b in bins[1:]:
    print b
    propweight = histdist[bincount]/np.sum(histdist)
    geoweight = 1./(b-binwidth/2.)**2.
    zmax = inv_comoving_distance(b)
    zmin = inv_comoving_distance(b-binwidth)
    slicearr = twomrsarr.transpose()[(z<zmax)*(z>zmin)].transpose()
    slicemap = scattomap(slicearr[1], slicearr[0])
    smoothangle = ((b-binwidth/2.)/10.)**0.5*100.*0.025*5./0.8 + 2.
    print 'Smoothing by', smoothangle
    slicemap = hp.smoothing(slicemap, sigma = np.deg2rad(smoothangle))
    slicemap = slicemap-1.*np.min(slicemap)
    slicemap = slicemap/slicemap.sum()
    avmap = avmap+propweight*slicemap
    
    bincount = bincount+1

avmap = hp.smoothing(avmap, sigma = np.deg2rad(4.))
avmap = avmap-1.*np.min(avmap)
avmap = avmap/avmap.sum()

hp.mollview(avmap, title="Weighted 2MRS. Heavier Composition")
plt.savefig('2MRSProjectedmapN__0.3N_0.3C_0.25Li_0.1He_0.05H_Be-9.png')
plt.show()
np.savetxt('2MRSProjectedmapN__0.3N_0.3C_0.25Li_0.1He_0.05H_Be-9.txt', avmap, delimiter="|")

#plt.hist(np.asarray(lengths))

#plt.show()




