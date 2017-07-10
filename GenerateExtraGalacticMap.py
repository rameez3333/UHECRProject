import matplotlib
matplotlib.use('Agg')
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import itertools
from scipy.interpolate import InterpolatedUnivariateSpline
from astropy.cosmology import Planck15 as cosmo
import math
from optparse import OptionParser
from wquantiles import quantile
from astropy import units as u
from astropy.coordinates import SkyCoord
#import multiprocessing
#pool = multiprocessing.Pool()

NSIDE=32

usage = 'usage: %prog [options]'
parser = OptionParser(usage)

parser.add_option("-a", "--massnumber", action="store", type="int", default=56, dest="MASSNUMBER", help="The Atomic Mass Number (A) of the Nucleus")
parser.add_option("-z", "--atomicnumber", action="store", type="int", default=26, dest="ATOMICNUMBER", help="The Atomic Number (Z) of the Nucleus")
parser.add_option("-l", "--zslices", action="store", type="int", default=30, dest="ZSLICES", help="The number of Z slices")
parser.add_option("-n", "--nsamples", action="store", type="int", default=1, dest="NSAMPLES", help="The number of MC maps to make")
parser.add_option("-i", "--ioffset", action="store", type="int", default=0, dest="IOFFSET", help="Index Offset")
parser.add_option("-d", "--density", action="store", type="float", default=1.e-4, dest="DENSITYSOUCE", help="UHECR source density in/MPc^3")
parser.add_option("-m", "--evolution", action="store", type="float", default=3.0, dest="MEVOL", help="Index of UHECR evolution")
parser.add_option("-b", "--magneticfield", action="store", type="float", default=1.e-10, dest="MAGFIELD", help="The extragalactic magnetic field (in Gauss)")
parser.add_option("-r", "--radius", action="store", type="float", default=300., dest="RADIUS", help="Radius of the sphercal volume to consider (MPc)")
parser.add_option("-s", "--spectralindex", action="store", type="float", default=2.0, dest="SPECTRALINDEX", help="The Spectral Index at Injection")
parser.add_option("-p", "--plots", action = "store_true", default=False, dest="PLOTS", help = "Whether to make plots")
parser.add_option("-o", "--otpf", action = "store_true", default=True, dest="OTPF", help = "One time plots and figures?")
parser.add_option("-f", "--lumfunc", action = "store_true", default=False, dest="LUMFUNC", help = "Use a luminosity function instead of standard candle sources?")
(options, args) = parser.parse_args()


comoving_distance = InterpolatedUnivariateSpline(np.linspace(0, 0.1, 100), cosmo.comoving_distance(np.linspace(0, 0.1, 100)).value)
inv_comoving_distance = InterpolatedUnivariateSpline(cosmo.comoving_distance(np.linspace(0, 0.1, 100)).value, np.linspace(0, 0.1, 100))

colours = ['red', 'blue', 'green', 'yellow', 'brown', 'purple', 'magenta', 'cyan', 'black']

ethresh = 53.0

galpixindexlist=None
restindexlist=None

def scattomap(dec,ra, nside=NSIDE, weights=None):
    hmap = np.bincount(hp.ang2pix(nside, np.deg2rad(90.-dec), np.deg2rad(ra)), minlength=hp.nside2npix(nside), weights=weights)
    return hmap

A = options.MASSNUMBER
Z = options.ATOMICNUMBER
Density = options.DENSITYSOUCE
Radius = options.RADIUS
mevol = options.MEVOL
zslices = options.ZSLICES
spectralindex = options.SPECTRALINDEX
magfield = options.MAGFIELD
nsamples=options.NSAMPLES
plots=options.PLOTS
otpf = options.OTPF
offset = options.IOFFSET
lumfunc = options.LUMFUNC

fnamebase = 'ProjectedMaps6/EGal_Injection_'+str(A)+'_'+str(Z)+'_SpecIndex_'+str(spectralindex)+'_Den_'+str(Density)+'_EgalMag_'+str(magfield)+'_'

if lumfunc:
    fnamebase = fnamebase+'LumFunc_'

Volume = 4./3.*np.pi*np.power(Radius, 3.)
nsources = int(Density*Volume)
zmax = inv_comoving_distance(Radius)
coherencelength = 10. #in Mpc

def IndexToDeclRa(index):
    theta,phi=hp.pixelfunc.pix2ang(NSIDE,index)
    return -np.degrees(theta-np.pi/2.),np.degrees(phi)

def AZtoNCode(A, Z):
    return 1e9 + Z*1e4 + A*1e1

def CodeToAZ(code):
    return int(str(code)[4:6]), int(str(code)[7:9])

class TwoMRSMap:
    def __init__(self, rawmap, mtype="SourceDist", smoothing = 0.0):
        global galpixindexlist
        global restindexlist
        if  galpixindexlist==None:
            indices = np.arange(hp.nside2npix(NSIDE))
            theta, phi = hp.pix2ang(NSIDE, indices)
            dec, ra = 90. - np.rad2deg(theta), np.rad2deg(phi)
            galb = SkyCoord(ra = ra*u.degree, dec=dec*u.degree).galactic.b.value
            galpixindexlist = indices[np.absolute(galb)<5.]
            restindexlist = indices[np.absolute(galb)>=5.]
            print 'Total to Masked ratio', len(indices) / len(galpixindexlist)
        if mtype=="SourceDist":
            self.EqNormMap, self.PosMap=self.MapToEqNorm(rawmap, smoothing)
        elif mtype=="TwoMRSDist":
            rawmap[galpixindexlist] = rawmap[np.random.choice(restindexlist, len(galpixindexlist))]
            self.EqNormMap, self.PosMap=self.MapToEqNorm(rawmap, smoothing)
        else:
            print PF("err"), "mtype==\""+mtype+"\" unknown"
            exit(1)
        self.SigWeight={}


    def MapToEqNorm(self,eqmap, smoothing):
        print 'wtf', eqmap, np.sum(eqmap) 
        if NSIDE!=32:
            if DEBUG:
                print PF("warn"),"rebinning EqNormMap to ",NSIDE
            eqmap=hp.pixelfunc.ud_grade(eqmap,NSIDE)
        
        if not np.sum(eqmap):
            return np.zeros(len(eqmap))
        if smoothing:
            eqmap = hp.smoothing(eqmap, sigma=np.deg2rad(smoothing),verbose=False,lmax=64)
        
        posmap = eqmap-1.*np.min(eqmap)
        eqmap[eqmap<0.]=0.
        if not np.sum(eqmap):
            eqmap = np.ones(len(eqmap))
        return eqmap/float(np.sum(eqmap)), posmap/float(np.sum(posmap))


def Make2MRSMap(redshiftcutlow=0.00, redshiftcuthigh=0.173, nside=NSIDE):
    twomrsarr = np.genfromtxt('2MRS/catalog/2mrs_1175_done.dat', skip_header=10, usecols=(1,2,24)).transpose()
    twomrsarr = np.vstack((twomrsarr[0], twomrsarr[1], twomrsarr[2]/299792.458))
    twomrsarr = twomrsarr.transpose()[twomrsarr[2]<redshiftcuthigh].transpose()
    twomrsarr = twomrsarr.transpose()[twomrsarr[2]>redshiftcutlow].transpose()
    return scattomap(twomrsarr[1], twomrsarr[0], nside)


print 'Generating ', nsources, ' sources with evolution m= ', mevol

Evoprobx = np.linspace(0, zmax, 1500)
Evoproby = np.power(1.+Evoprobx, mevol)*np.power(Evoprobx, 2.)
backpolate = InterpolatedUnivariateSpline(Evoproby, Evoprobx)


Zarr = np.linspace(0., zmax, zslices)


TomoMaps={}
for i in range(len(Zarr)-1):
    print 'Loading map in range ', Zarr[i], ' to ', Zarr[i+1]
    TomoMaps[float('%.4g' % Zarr[i])] = TwoMRSMap(Make2MRSMap(Zarr[i],  Zarr[i+1], NSIDE), "TwoMRSDist", 0.)

#print sorted(TomoMaps.keys())

def GenerateDecRaReal(z):
        if z<Zarr[-2]:
            znearest = Zarr[np.abs(Zarr - z).argmin()]
            znearest = float('%.4g' % znearest)
            samplemap = TomoMaps[znearest].EqNormMap
        else:
            samplemap = np.ones(hp.nside2npix(NSIDE))/float(hp.nside2npix(NSIDE))
        
        index = np.random.choice(np.arange(hp.nside2npix(NSIDE)), size=1, p = samplemap)
        
        return IndexToDeclRa(index[0])
    
def GeneRateSources(nsources):
    rs = np.random.uniform(0, np.power((1.+zmax), mevol)*np.power(zmax,2), size=nsources)
    redshifts = backpolate(rs)
    decra = np.asarray(map(GenerateDecRaReal, redshifts)).transpose()
    distances = comoving_distance(redshifts)
    decs, ras =  decra[0], decra[1]
    return np.vstack((ras, decs, redshifts))


def MakeGeneratedMap(twomrsarr, redshiftcutlow=0.00, redshiftcuthigh=0.173, nside=NSIDE, lumfunc=False):
    #twomrsarr = np.vstack((ras, decs, redshifts))
    twomrsarr = twomrsarr.transpose()[twomrsarr[2]<redshiftcuthigh].transpose()
    twomrsarr = twomrsarr.transpose()[twomrsarr[2]>redshiftcutlow].transpose()
    weights=None
    if lumfunc:
        weights = np.power(10., np.random.normal(0.0, 0.4472, len(twomrsarr[1])))
    return scattomap(twomrsarr[1], twomrsarr[0], nside, weights=weights)

#print len(decs), len(ras), nsources

#twomrsarr = np.genfromtxt('2MRS/catalog/2mrs_1175_done.dat', skip_header=10, usecols=(1,2,24)).transpose()
#twomrstype = np.genfromtxt('2MRS/catalog/2mrs_1175_done.dat', skip_header=10, usecols=(22), dtype="|S5").transpose()
#twomrsmag = np.genfromtxt('2MRS/catalog/2mrs_1175_done.dat', skip_header=10, usecols=(8,9,10)).transpose()
#print len(twomrsarr[0]), 'objects loaded from the 2MRS catalog'
#z = twomrsarr[2]/299792.458




#plt.hist(twomrsmag[0])
#plt.show()

#print 'Selecting Magnitude < 10'
#twomrsarr = twomrsarr.transpose()[twomrsmag[0]<10.].transpose()

#twomrsarr = twomrsarr.transpose()[z<0.01].transpose()
#z = twomrsarr[2]/299792.458




#del plt
#import matplotlib.pyplot as plt

#plt.figure(2)
#map2 = scattomap(twomrsarr[1], twomrsarr[0])
#hp.mollview(map2, title='2MRS Sources')
#plt.savefig('2MRS.png')


print 'Loading File'

propagated = np.genfromtxt('Outputs/'+str(A)+'_'+str(Z)+'/MinEn53.0MaxEn400.0_InSpecIndex'+str(spectralindex)+'_StartDist500.0_N1000000_Seq0.txt')


print 'Loaded'
startindices = np.where(propagated.transpose()[0]==0.)[0]
weights = np.insert(startindices, len(startindices), len(propagated))
weights  = weights[1:] - startindices
counts = weights
weights = weights*1.0
weights = 1./weights

print 'maps'
orig = np.asarray(map(CodeToAZ, propagated.transpose()[3]))
det = np.asarray(map(CodeToAZ, propagated.transpose()[1]))
print 'mapped'

origz, origa = orig.transpose()[1], orig.transpose()[0]
detz, deta = det.transpose()[1], det.transpose()[0]

propagated = propagated.transpose()

propagated = np.vstack([propagated[0], deta, detz, propagated[2], propagated[4]]).transpose()
weights = np.repeat(weights, counts)

weights = weights[propagated.transpose()[3]>ethresh]
propagated=propagated[propagated.transpose()[3]>ethresh]

#histaz = np.histogram2d(propagated.transpose()[1], propagated.transpose()[2], weights=weights, bins=[27,56], range=[[0,27],[1,57]])
#zbins = histaz[1]
#abins = histaz[2]
#histazobtained = histaz[0]/np.sum(histaz[0])

#distatthreshold={}

#zindices = np.digitize(propagated.transpose()[1], zbins, right=True)
#aindices = np.digitize(propagated.transpose()[2], abins, right=True)

#distatinjection={}

#stepsize = 2

#detectedenergies = np.arange(54, 200, stepsize)
#injectedenergies = np.arange(58, 198, stepsize)

#meandistatdetection=[]
#meandistatinjection=[]
weighttot = np.sum(weights)

for count in range(0+offset,nsamples+offset):
    fname=fnamebase+'S'+str(count)+'_'
    
    print 'F: ', fnamebase
    
    Nsourceslice={}
    histazslice={}
    GenerateMaps={}
    arrmap = np.zeros(hp.nside2npix(NSIDE))
    meddistmap = np.zeros(hp.nside2npix(NSIDE))
    totmap = np.zeros(hp.nside2npix(NSIDE))
    twomrsarr = GeneRateSources(nsources)
    comphistlist=[]
    medano=[]
    ano1=[]
    ano5=[]
    ano9=[]
    meden=[]
    en1=[]
    en5=[]
    en9=[]
    if plots:
        print 'Making Plots'
        ras, decs = twomrsarr[0], twomrsarr[1]
        plt.figure(1)
        map1 = scattomap(decs, ras)
        hp.mollview(map1, title='Generated Sources')
        plt.savefig(fname+'GeneratedSources.png')

    redshifts = twomrsarr[2]
    for i in range(len(Zarr)-1):
        #try:
        print 'Now processing ', Zarr[i], ' to ', Zarr[i+1], 'i.e', comoving_distance(Zarr[i]), ' to ', comoving_distance(Zarr[i+1])
        meddist = (comoving_distance(Zarr[i]) + comoving_distance(Zarr[i+1]))/2.
        Nsourceslice = len(redshifts[(redshifts>Zarr[i])*(redshifts<Zarr[i+1])])
        print 'Sources in slice, total ',   Nsourceslice, nsources
        slicepropselect = propagated[(propagated.transpose()[0]>comoving_distance(Zarr[i]))*(propagated.transpose()[0]<comoving_distance(Zarr[i+1]))]
        sliceweights = weights[(propagated.transpose()[0]>comoving_distance(Zarr[i]))*(propagated.transpose()[0]<comoving_distance(Zarr[i+1]))]
        histazslice = np.histogram2d(slicepropselect.transpose()[1], slicepropselect.transpose()[2], weights=sliceweights, bins=[27,56], range=[[0,27],[1,57]])
        histazslice = histazslice[0]/np.sum(histazslice[0])
        comphistlist.append(histazslice)
        medenergy = np.sum(slicepropselect.transpose()[3]*sliceweights)/np.sum(sliceweights)
        energy1=quantile(slicepropselect.transpose()[3], sliceweights, 0.1)
        energy5=quantile(slicepropselect.transpose()[3], sliceweights, 0.5)
        energy9=quantile(slicepropselect.transpose()[3], sliceweights, 0.9)
        medatonum = np.sum(slicepropselect.transpose()[1]*sliceweights)/np.sum(sliceweights)
        atonu1=quantile(slicepropselect.transpose()[1], sliceweights, 0.1)
        atonu5=quantile(slicepropselect.transpose()[1], sliceweights, 0.5)
        atonu9=quantile(slicepropselect.transpose()[1], sliceweights, 0.9)
        print 'Mean Energy from this slice', medenergy
        #print 'Composition from this slice', histazslice
        print 'Mean atomic number from this slice', medatonum
        medano.append(medatonum)
        ano1.append(atonu1)
        ano5.append(atonu5)
        ano9.append(atonu9)
        meden.append(medenergy)
        en1.append(energy1)
        en5.append(energy5)
        en9.append(energy9)
        Zato = (medatonum + Z)/2.
        if not Nsourceslice:
            continue
        defang = 0.025*np.sqrt(meddist/coherencelength)*(coherencelength/10.)*(magfield/1e-11)*(100./medenergy)*Zato
        print 'Deflection ', defang    
        GenerateMaps = TwoMRSMap(MakeGeneratedMap(twomrsarr, Zarr[i],  Zarr[i+1], NSIDE, lumfunc), "SourceDist", defang)
        weight = float(Nsourceslice)/float(nsources)/4./np.pi/meddist**2.
        extweight = np.sum(sliceweights)/weighttot
        print 'weight due to extinction', extweight
        weight = weight*extweight
        arrmap = arrmap + weight*GenerateMaps.EqNormMap
        totmap = totmap + GenerateMaps.PosMap
        if np.any(np.isnan(GenerateMaps.EqNormMap)):
            print 'Nans in EqNormMap'
        if np.any(np.isnan(GenerateMaps.EqNormMap)):
            print 'Nans in PosMap'
        meddistmap = meddistmap + meddist*GenerateMaps.PosMap
        print meddist, weight,  GenerateMaps.EqNormMap
        #except:
            #print 'Stopping at slice no :', i
            #break
        
    if np.any(np.isnan(1./totmap)):
        print 'NaNs in Totmap inverse'
        print np.min(totmap)
    arrmap = arrmap/np.sum(arrmap)
    meddistmap = meddistmap/totmap
    
    if np.any(np.isnan(meddistmap)):
        meddistmap[np.isnan(meddistmap)] = Radius
    
    savemap = np.vstack([arrmap, meddistmap])
    if np.any(np.isnan(np.vstack([arrmap, meddistmap]))):
        print 'NaNs in Map'
    np.savetxt(fname+'Map.txt', savemap, delimiter="|")
    
    if (otpf) and (not count):
        np.savetxt(fname+'CompHist.txt', np.vstack(comphistlist), delimiter="|")
        print len(Zarr), len(medano), len(ano1), len(ano5), len(ano9), len(meden), len(en1), len(en5), len(en9), i
        np.savetxt(fname+'CompEn.txt', np.vstack([Zarr[0:i-1], medano[0:i-1], ano1[0:i-1], ano5[0:i-1], ano9[0:i-1], meden[0:i-1], en1[0:i-1], en5[0:i-1], en9[0:i-1]]), delimiter="|")
        if np.any(np.isnan(np.vstack([Zarr[0:i-1], medano[0:i-1], ano1[0:i-1], ano5[0:i-1], ano9[0:i-1], meden[0:i-1], en1[0:i-1], en5[0:i-1], en9[0:i-1]]))):
            print 'NaNs in CompEn'
    
    
    
    if plots:
        plt.figure(2)
        hp.mollview(arrmap)
        plt.savefig(fname+'Map.png')
        plt.figure(3)
        hp.mollview(meddistmap)
        plt.savefig(fname+'meddistMap.png')
    #plt.show()
    
    
    #TomoMaps[float('%.4g' % Zarr[i])]



#for denergy in detectedenergies:
    #print denergy
    #distatthreshold[denergy]=[]
    #for i in xrange(0,len(startindices)-1):
        #distatthreshold[denergy].append(propagated[startindices[i]:startindices[i+1]].transpose()[0][np.where(propagated[startindices[i]:startindices[i+1]].transpose()[2]>denergy)][-1])
    #meandistatdetection.append(np.median(distatthreshold[denergy]))

#for ienergy in injectedenergies:
    #print ienergy
    #distatinjection[ienergy]=[]
    #for i in xrange(0,len(startindices)-1):
        #distatinjection[ienergy].append(propagated[startindices[i]:startindices[i+1]].transpose()[0][np.where(propagated[startindices[i]:startindices[i+1]].transpose()[2]>ienergy)][-1])
    #meandistatinjection.append(np.mean(distatinjection[ienergy]))
    
    

#meandistatdetection is for a proton injected at 200 EeV

#dist = InterpolatedUnivariateSpline(detectedenergies, meandistatdetection)
    
##plt.plot(injectedenergies, meandistatinjection)

#plt.plot(detectedenergies, meandistatdetection, color='green', linewidth=2, alpha=0.5)
#plt.show()

#def generatepowerlaw(index, rmin, rmax, size):
    #arr = np.random.uniform(np.power(rmin, -1.*index), np.power(rmax, -1.*index), size)
    #return np.power(arr, -1./index)
    

#detectedspectralindex = 4.3

#injectionspectralindices = np.arange(0.8, 2.2, 0.2)
#energiesatdetection = generatepowerlaw(detectedspectralindex, 57., 190., 1000)
#sampleddistances = {}
#i=0
#for injectionspectralindex in injectionspectralindices:
    #energiesatinjection = generatepowerlaw(injectionspectralindex, 57., 200., 1000)


    #particles = [(x, y) for x in energiesatdetection for y in energiesatinjection]

    ##distances=[]


    #denbins = detectedenergies[map(lambda x : (np.abs(detectedenergies - x[0])).argmin(), particles)]
    #ienbins = detectedenergies[map(lambda x : (np.abs(detectedenergies - x[1])).argmin(), particles)]

    ##bins = np.arange(50, 200, 10)

    ##plt.hist(denbins, color='blue', alpha=0.5, bins=bins)
    ##plt.hist(ienbins, color='blue', alpha=0.5, bins=bins)
    ##plt.xscale('log')
    ##plt.yscale('log')
    ##plt.show()

    #di = np.random.randint(0, len(distatthreshold[104]), size = len(particles))
    #ii = np.random.randint(0, len(distatthreshold[104]), size = len(particles))

    #dsample = np.asarray(map(lambda p,q : distatthreshold[p][q], denbins, di))
    #isample = np.asarray(map(lambda p,q : distatthreshold[p][q], ienbins, ii))

    #distances = dsample - isample

    #sampleddistances[injectionspectralindex] = distances[distances>0]

    ##for den, ien in particles:
        ##if den <= ien:
            ##print den, ien
            ##distances.append(dist(den) - dist(ien))
            ##denbin = find_nearest(detectedenergies, den)
            ##ienbin = find_nearest(detectedenergies, ien)
            ##dsample = np.random.choice(distatthreshold[denbin])
            ##isample = np.random.choice(distatthreshold[ienbin])
            ##if (dsample - isample) > 0:
                ##sampleddistances.append(dsample-isample)

    #medredshift = inv_comoving_distance(np.median(sampleddistances[injectionspectralindex]))
            
    ##print 'Median Distance', np.median(distances)
    #print 'Spectral Index', injectionspectralindex
    #print 'Median distance', np.median(sampleddistances[injectionspectralindex])
    #print 'Median Redshift', inv_comoving_distance(np.median(sampleddistances[injectionspectralindex]))

    #print len(twomrsarr.transpose()[z<medredshift].transpose()[0]), '2MRS sources within this distance'

    #bins = np.arange(0, 600, 2)

    ##plt.hist(distances, color='blue', alpha=0.5, bins=bins)
    #plt.hist(sampleddistances[injectionspectralindex], color=colours[i], alpha=0.2, bins=bins, label=str(injectionspectralindex))
    #i=i+1
#plt.legend(loc='best', fontsize=15) 
#plt.xlabel('Distance (MPc)')
#plt.ylabel('Events')
#plt.yscale('log')
    ##plt.show()

#plt.savefig('DistrosProtonOnly.png')
        
#plt.show()

#del plt

#import matplotlib.pyplot as plt

#def scattomap(dec,ra, nside=32):
    #hmap = np.bincount(hp.ang2pix(nside, np.deg2rad(90.-dec), np.deg2rad(ra)), minlength=hp.nside2npix(nside))
    #return hmap


#injectionspectralindex = injectionspectralindices[-1]
#print 'Injection Index', injectionspectralindex

#avmap = np.zeros(hp.nside2npix(32))

#ncount=[]
#for b in bins[1:]:
    #print b
    #propweight = float(len(sampleddistances[injectionspectralindex][(sampleddistances[injectionspectralindex] < b)*(sampleddistances[injectionspectralindex] > (b-10.))]))/float(len(sampleddistances[injectionspectralindex]))
    #geoweight = 1./(b-1.)**2.
    #zmax = inv_comoving_distance(b)
    #zmin = inv_comoving_distance(b-2.)
    #slicearr = twomrsarr.transpose()[(z<zmax)*(z>zmin)].transpose()
    #ncount.append(len(slicearr[0]))
    #slicemap = scattomap(slicearr[1], slicearr[0])
    ##hp.mollview(slicemap)
    ##plt.show()
    #print slicemap
    #print avmap
    #print propweight
    #hp.mollview(slicemap, title=str(b)+ ' MPc')
    #plt.savefig('2MRSMaps/map'+str(len(ncount))+'.png')
    #if b<120:
        #avmap = avmap+propweight*geoweight*slicemap
    #elif ncount[-1]:
        #avmap = avmap+propweight*geoweight*slicemap*float(ncount[20])/float(ncount[-1])
    #else:
        #continue
    

#avmap = avmap/np.sum(avmap)

#plt.plot(bins[1:], ncount)
#plt.xlabel('Distance (MPc)')
#plt.ylabel('NSources in 2MRS within 2MPc shells')
#plt.show()

#hp.mollview(avmap, title="Weighted 2MRS")
#plt.savefig('2MRSProjectedMapPureProtonFineBin.png')
#plt.show()
#np.savetxt('2MRSProjectedMapPureProtonFineBin.txt', avmap, delimiter="|")
