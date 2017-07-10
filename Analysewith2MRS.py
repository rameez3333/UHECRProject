import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import itertools
from scipy.interpolate import InterpolatedUnivariateSpline
from astropy.cosmology import Planck15 as cosmo
import math


inv_comoving_distance = InterpolatedUnivariateSpline(cosmo.comoving_distance(np.linspace(0, 0.1, 100)).value, np.linspace(0, 0.1, 100))



colours = ['red', 'blue', 'green', 'yellow', 'brown', 'purple', 'magenta', 'cyan', 'black']

def scattomap(dec,ra, nside=32):
    hmap = np.bincount(hp.ang2pix(nside, np.deg2rad(90.-dec), np.deg2rad(ra)), minlength=hp.nside2npix(nside))
    return hmap

A = 1
Z = 1

twomrsarr = np.genfromtxt('2MRS/catalog/2mrs_1175_done.dat', skip_header=10, usecols=(1,2,24)).transpose()
twomrstype = np.genfromtxt('2MRS/catalog/2mrs_1175_done.dat', skip_header=10, usecols=(22), dtype="|S5").transpose()
twomrsmag = np.genfromtxt('2MRS/catalog/2mrs_1175_done.dat', skip_header=10, usecols=(8,9,10)).transpose()
print len(twomrsarr[0]), 'objects loaded from the 2MRS catalog'
z = twomrsarr[2]/299792.458

#plt.hist(twomrsmag[0])
#plt.show()

#print 'Selecting Magnitude < 10'
#twomrsarr = twomrsarr.transpose()[twomrsmag[0]<10.].transpose()

#twomrsarr = twomrsarr.transpose()[z<0.01].transpose()
z = twomrsarr[2]/299792.458

map1 = scattomap(twomrsarr[1], twomrsarr[0])
#hp.mollview(map1)
#plt.show()

#print 'Selecting AGNs'
#y = np.asarray(map(lambda x : True if x[0]=='-' else False, twomrstype))
#twomrsarr = twomrsarr.transpose()[y].transpose()
#twomrstype = twomrstype[y]
#y = np.asarray(map(lambda x : True if x[1]=='9' else False, twomrstype))
#twomrsarr = twomrsarr.transpose()[y].transpose()
#twomrstype = twomrstype[y]
#z = twomrsarr[2]/299792.458

propagated = np.genfromtxt('Outputs/'+str(A)+'_'+str(Z)+'/MinEn53.0_StartEn200.0_StartDist200.0_N1000000_Seq0.txt')

startindices = np.where(propagated.transpose()[0]==0.)[0]

distatthreshold={}
#distatinjection={}

stepsize = 2

detectedenergies = np.arange(54, 200, stepsize)
#injectedenergies = np.arange(58, 198, stepsize)

meandistatdetection=[]
#meandistatinjection=[]


def AZtoNCode(A, Z):
    return 1e9 + Z*1e4 + A*1e1


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

#plt.plot(detectedenergies, meandistatdetection, color='green', linewidth=2, alpha=0.5)
#plt.show()

def generatepowerlaw(index, rmin, rmax, size):
    arr = np.random.uniform(np.power(rmin, -1.*index), np.power(rmax, -1.*index), size)
    return np.power(arr, -1./index)
    

detectedspectralindex = 4.3

injectionspectralindices = np.arange(0.8, 2.2, 0.2)
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

    bins = np.arange(0, 600, 2)

    #plt.hist(distances, color='blue', alpha=0.5, bins=bins)
    #plt.hist(sampleddistances[injectionspectralindex], color=colours[i], alpha=0.2, bins=bins, label=str(injectionspectralindex))
    i=i+1
#plt.legend(loc='best', fontsize=15) 
#plt.xlabel('Distance (Mpc)')
#plt.ylabel('Events')
#plt.yscale('log')
    #plt.show()

#plt.savefig('DistrosProtonOnly.png')
        
#plt.show()

del plt

import matplotlib.pyplot as plt

def scattomap(dec,ra, nside=32):
    hmap = np.bincount(hp.ang2pix(nside, np.deg2rad(90.-dec), np.deg2rad(ra)), minlength=hp.nside2npix(nside))
    return hmap


injectionspectralindex = injectionspectralindices[-1]
print 'Injection Index', injectionspectralindex

avmap = np.zeros(hp.nside2npix(32))

ncount=[]
for b in bins[1:]:
    print b
    propweight = float(len(sampleddistances[injectionspectralindex][(sampleddistances[injectionspectralindex] < b)*(sampleddistances[injectionspectralindex] > (b-10.))]))/float(len(sampleddistances[injectionspectralindex]))
    geoweight = 1./(b-1.)**2.
    zmax = inv_comoving_distance(b)
    zmin = inv_comoving_distance(b-2.)
    slicearr = twomrsarr.transpose()[(z<zmax)*(z>zmin)].transpose()
    ncount.append(len(slicearr[0]))
    slicemap = scattomap(slicearr[1], slicearr[0])
    #hp.mollview(slicemap)
    #plt.show()
    print slicemap
    print avmap
    print propweight
    #hp.mollview(slicemap, title=str(b)+ ' MPc')
    plt.savefig('2MRSMaps/map'+str(len(ncount))+'.png')
    if b<120:
        avmap = avmap+propweight*geoweight*slicemap
    elif ncount[-1]:
        avmap = avmap+propweight*geoweight*slicemap*float(ncount[20])/float(ncount[-1])
    else:
        continue
    

avmap = avmap/np.sum(avmap)

plt.plot(bins[1:], ncount)
plt.xlabel('Distance (Mpc)', fontsize=15)
plt.ylabel('No of Sources in 2MRS within 2Mpc shells', fontsize=15)
plt.savefig('NsourcesVsDist.png')
plt.show()

#hp.mollview(avmap, title="Weighted 2MRS")
#plt.savefig('2MRSProjectedMapPureProtonFineBin.png')
#plt.show()
#np.savetxt('2MRSProjectedMapPureProtonFineBin.txt', avmap, delimiter="|")
