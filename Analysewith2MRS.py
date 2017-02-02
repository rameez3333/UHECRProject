import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import itertools
from scipy.interpolate import InterpolatedUnivariateSpline
from astropy.cosmology import Planck15 as cosmo
import math


inv_comoving_distance = InterpolatedUnivariateSpline(cosmo.comoving_distance(np.linspace(0, 0.1, 100)).value, np.linspace(0, 0.1, 100))



colours = ['red', 'blue', 'green', 'yellow', 'brown', 'purple']


twomrsarr = np.genfromtxt('2MRS/catalog/2mrs_1175_done.dat', skip_header=10, usecols=(1,2,24)).transpose()
print len(twomrsarr[0]), 'objects loaded from the 2MRS catalog'
z = twomrsarr[2]/299792.458



propagated = np.genfromtxt('Outputs/1_1/MinEn53.0_StartEn200.0_StartDist100.0_N100000_Seq2.txt')

startindices = np.where(propagated.transpose()[0]==0.)[0]

distatthreshold={}
#distatinjection={}

stepsize = 2

detectedenergies = np.arange(54, 200, stepsize)
#injectedenergies = np.arange(58, 198, stepsize)

meandistatdetection=[]
#meandistatinjection=[]


for denergy in detectedenergies:
    print denergy
    distatthreshold[denergy]=[]
    for i in xrange(0,len(startindices)-1):
        distatthreshold[denergy].append(propagated[startindices[i]:startindices[i+1]].transpose()[0][np.where(propagated[startindices[i]:startindices[i+1]].transpose()[2]>denergy)][-1])
    meandistatdetection.append(np.median(distatthreshold[denergy]))

#for ienergy in injectedenergies:
    #print ienergy
    #distatinjection[ienergy]=[]
    #for i in xrange(0,len(startindices)-1):
        #distatinjection[ienergy].append(propagated[startindices[i]:startindices[i+1]].transpose()[0][np.where(propagated[startindices[i]:startindices[i+1]].transpose()[2]>ienergy)][-1])
    #meandistatinjection.append(np.mean(distatinjection[ienergy]))
    
    

#meandistatdetection is for a proton injected at 200 EeV

dist = InterpolatedUnivariateSpline(detectedenergies, meandistatdetection)
    
#plt.plot(injectedenergies, meandistatinjection)

plt.plot(detectedenergies, meandistatdetection, color='green', linewidth=2, alpha=0.5)
plt.show()

def generatepowerlaw(index, rmin, rmax, size):
    arr = np.random.uniform(np.power(rmin, -1.*index), np.power(rmax, -1.*index), size)
    return np.power(arr, -1./index)
    

detectedspectralindex = 4.3

injectionspectralindices = np.arange(2., 3.1, 0.2)
energiesatdetection = generatepowerlaw(detectedspectralindex, 57., 190., 1000)
sampleddistances = {}
i=0
for injectionspectralindex in injectionspectralindices:
    energiesatinjection = generatepowerlaw(injectionspectralindex, 57., 200., 1000)


    particles = [(x, y) for x in energiesatdetection for y in energiesatinjection]

    #distances=[]


    denbins = detectedenergies[map(lambda x : (np.abs(detectedenergies - x[0])).argmin(), particles)]
    ienbins = detectedenergies[map(lambda x : (np.abs(detectedenergies - x[1])).argmin(), particles)]

    #bins = np.arange(50, 200, 10)

    #plt.hist(denbins, color='blue', alpha=0.5, bins=bins)
    #plt.hist(ienbins, color='blue', alpha=0.5, bins=bins)
    #plt.xscale('log')
    #plt.yscale('log')
    #plt.show()

    di = np.random.randint(0, len(distatthreshold[104]), size = len(particles))
    ii = np.random.randint(0, len(distatthreshold[104]), size = len(particles))

    dsample = np.asarray(map(lambda p,q : distatthreshold[p][q], denbins, di))
    isample = np.asarray(map(lambda p,q : distatthreshold[p][q], ienbins, ii))

    distances = dsample - isample

    sampleddistances[injectionspectralindex] = distances[distances>0]

    #for den, ien in particles:
        #if den <= ien:
            #print den, ien
            #distances.append(dist(den) - dist(ien))
            #denbin = find_nearest(detectedenergies, den)
            #ienbin = find_nearest(detectedenergies, ien)
            #dsample = np.random.choice(distatthreshold[denbin])
            #isample = np.random.choice(distatthreshold[ienbin])
            #if (dsample - isample) > 0:
                #sampleddistances.append(dsample-isample)

    medredshift = inv_comoving_distance(np.median(sampleddistances[injectionspectralindex]))
            
    #print 'Median Distance', np.median(distances)
    print 'Spectral Index', injectionspectralindex
    print 'Median distance', np.median(sampleddistances[injectionspectralindex])
    print 'Median Redshift', inv_comoving_distance(np.median(sampleddistances[injectionspectralindex]))

    print len(twomrsarr.transpose()[z<medredshift].transpose()[0]), '2MRS sources within this distance'

    bins = np.arange(0, 600, 10)

    #plt.hist(distances, color='blue', alpha=0.5, bins=bins)
    plt.hist(sampleddistances[injectionspectralindex], color=colours[i], alpha=0.2, bins=bins, label=str(injectionspectralindex))
    i=i+1
plt.legend(loc='best', fontsize=15)
plt.xlabel('Distance (MPc)')
plt.ylabel('Events')
plt.yscale('log')
plt.show()
plt.savefig('Distros.png')
        

def scattomap(dec,ra, nside=32):
    hmap = np.bincount(hp.ang2pix(nside, np.deg2rad(90.-dec), np.deg2rad(ra)), minlength=hp.nside2npix(nside))
    return hmap


injectionspectralindex = 2.

avmap = np.zeros(hp.nside2npix(32))

for b in bins[1:]:
    print b
    propweight = float(len(sampleddistances[injectionspectralindex][(sampleddistances[injectionspectralindex] < b)*(sampleddistances[injectionspectralindex] > (b-10.))]))/float(len(sampleddistances[injectionspectralindex]))
    geoweight = 1./(b-5.)**2.
    zmax = inv_comoving_distance(b)
    zmin = inv_comoving_distance(b-10.)
    slicearr = twomrsarr.transpose()[(z<zmax)*(z>zmin)].transpose()
    
    slicemap = scattomap(slicearr[1], slicearr[0])
    #hp.mollview(slicemap)
    #plt.show()
    print slicemap
    print avmap
    print propweight
    avmap = avmap+propweight*geoweight*slicemap

hp.mollview(avmap, title="Weighted 2MRS")
plt.show()
np.savetxt('2MRSProjectedmap.txt', avmap, delimiter="|")
