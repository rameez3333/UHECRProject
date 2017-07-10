import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import healpy as hp
from glob import glob
from scipy.interpolate import InterpolatedUnivariateSpline
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.cosmology import Planck15 as cosmo
import Data
import CR_funcs as CR
import os, sys
from optparse import OptionParser

usage = 'usage: %prog [options]'
parser = OptionParser(usage)

parser.add_option("-a", "--massnumber", action="store", type="int", default=56, dest="MASSNUMBER", help="The Atomic Mass Number (A) of the Nucleus")
parser.add_option("-z", "--atomicnumber", action="store", type="int", default=26, dest="ATOMICNUMBER", help="The Atomic Number (Z) of the Nucleus")
parser.add_option("-s", "--spectralindex", action="store", type="float", default=2.0, dest="SPECTRALINDEX", help="The Spectral Index at Injection")
parser.add_option("-d", "--density", action="store", type="float", default=1.e-4, dest="DENSITYSOUCE", help="UHECR source density in/MPc^3")
parser.add_option("-b", "--magneticfield", action="store", type="float", default=1.e-10, dest="MAGFIELD", help="The extragalactic magnetic field (in Gauss)")
parser.add_option("-p", "--PT2011", action = "store_true", default=False, dest="PT2011", help = "Use the Pshirkov 2011 Field instead of the JF2012 default?")

(options, args) = parser.parse_args()

A = options.MASSNUMBER
Z = options.ATOMICNUMBER
Density = options.DENSITYSOUCE
spectralindex = options.SPECTRALINDEX
magfield = options.MAGFIELD

fnamebase = 'ProjectedMaps5/EGal_Injection_'+str(A)+'_'+str(Z)+'_SpecIndex_'+str(spectralindex)+'_Den_'+str(Density)+'_EgalMag_'+str(magfield)+'_S*_Map.txt'


allflist = sorted(glob(fnamebase))

E_min_Auger = 57
E_scale = 1.13
events, ff = Data.load_Auger_TA(E_min_Auger = E_min_Auger, E_scale = E_scale)

nucs = [[1,1], [4, 2],[7, 3], [9,4], [11,5],[12, 6], [14,7],[16,8], [19, 9], [20, 10], [23, 11], [24, 12], [27,13], [28,14], [31, 15], [32, 16], [35, 17],[40,18],[39, 19], [40, 20],[45,21],[48,22],[51,23], [52, 24],[55,25], [56, 26]]

nucmassdict = {}

for nuc in nucs:
    nucmassdict[nuc[1]] = nuc[0]

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

def scattomap(dec,ra, nside=32):
    hmap = np.bincount(hp.ang2pix(nside, np.deg2rad(90.-dec), np.deg2rad(ra)), minlength=hp.nside2npix(nside))
    return hmap    

def weightedscattomap(dec,ra, weights, nside=32):
    hmap = np.bincount(hp.ang2pix(nside, np.deg2rad(90.-dec), np.deg2rad(ra)), weights=weights, minlength=hp.nside2npix(nside))
    return hmap

inv_comoving_distance = InterpolatedUnivariateSpline(cosmo.comoving_distance(np.linspace(0, 0.1, 100)).value, np.linspace(0, 0.1, 100))

Radius = 300.
zslices = 30.


zmax = inv_comoving_distance(Radius)
Zarr = np.linspace(0., zmax, zslices)


comps = [[1,1], [16,8], [28,14], [56,26]]
specs = [0.9, 2.0]
densities = [1.e-5, 1.e-4, 1.e-3]
magfields = [1.e-11, 1.e-10, 1.e-9]

otfarr={}
compmeddict={}
comp1dict={}
comp5dict={}
comp9dict={}
enmeddict={}
en1dict={}
en5dict={}
en9dict={}
figcount=0
for comp in comps:
    for spec in specs:
        for density in densities:
            for mag in magfields:
                try:
                    figcount=figcount+1
                    otfarr[str(comp[0])+'_'+str(spec)+'_'+str(density)+str(mag)] = np.genfromtxt('ProjectedMaps/EGal_Injection_'+str(comp[0])+'_'+str(comp[1])+'_SpecIndex_'+str(spec)+'_Den_' + str(density) + '_EgalMag_' +str(mag)+'_S0_CompEn.txt', delimiter="|")
                    arrnow = otfarr[str(comp[0])+'_'+str(spec)+'_'+str(density)+str(mag)]
                    zarrnow = Zarr[0:len(arrnow[0])]
                    zwidth = Zarr[1]-Zarr[0]
                    zarrnow = zarrnow+zwidth/2.
                    #plt.figure(figcount)
                    compmeddict[str(comp[0])+'_'+str(spec)+'_'+str(density)+str(mag)] = InterpolatedUnivariateSpline(zarrnow, arrnow[1])
                    comp1dict[str(comp[0])+'_'+str(spec)+'_'+str(density)+str(mag)] = InterpolatedUnivariateSpline(zarrnow, arrnow[2])
                    comp5dict[str(comp[0])+'_'+str(spec)+'_'+str(density)+str(mag)] = InterpolatedUnivariateSpline(zarrnow, arrnow[3])
                    comp9dict[str(comp[0])+'_'+str(spec)+'_'+str(density)+str(mag)] = InterpolatedUnivariateSpline(zarrnow, arrnow[4])
                    enmeddict[str(comp[0])+'_'+str(spec)+'_'+str(density)+str(mag)] = InterpolatedUnivariateSpline(zarrnow, arrnow[5])
                    en1dict[str(comp[0])+'_'+str(spec)+'_'+str(density)+str(mag)] = InterpolatedUnivariateSpline(zarrnow, arrnow[6])
                    en5dict[str(comp[0])+'_'+str(spec)+'_'+str(density)+str(mag)] = InterpolatedUnivariateSpline(zarrnow, arrnow[7])
                    en9dict[str(comp[0])+'_'+str(spec)+'_'+str(density)+str(mag)] = InterpolatedUnivariateSpline(zarrnow, arrnow[8])
                    
                    #plt.plot(zarrnow, arrnow[1], label='mean', color='black')
                    #plt.plot(zarrnow, arrnow[2], label='10', color='red')
                    #plt.plot(zarrnow, arrnow[3], label='50', color='green')
                    #plt.plot(zarrnow, arrnow[4], label='90', color='blue')
                    #plt.xlabel('redshift')
                    #plt.ylabel('Atomic Number')
                    #plt.legend(loc='best')
                    #plt.savefig('Z_S_D_B_'+str(comp[0])+'_'+str(spec)+'_'+str(density)+str(mag)+'Comp.png')
                    #figcount=figcount+1
                    #plt.figure(figcount)
                    #plt.plot(zarrnow, arrnow[5], label='mean', color='black')
                    #plt.plot(zarrnow, arrnow[6], label='10', color='red')
                    #plt.plot(zarrnow, arrnow[7], label='50', color='green')
                    #plt.plot(zarrnow, arrnow[8], label='90', color='blue')
                    #plt.xlabel('redshift')
                    #plt.ylabel('Energy [EeV]')
                    #plt.legend(loc='best')
                    #plt.savefig('Z_S_D_B_'+str(comp[0])+'_'+str(spec)+'_'+str(density)+str(mag)+'En.png')
                except:
                    continue
        

evdict={}
evweightdict={}

def ProcessMap(fname):
    print 'Processing :', fname
    import matplotlib.pyplot as plt
    A = int(fname[fname.find('ion_')+4:fname.find('ion_')+6].replace('_',''))
    Z = int(fname[fname.find(str(A))+2:fname.find(str(A))+len(str(A))+3].replace('_',''))
    spec = float(fname[fname.find('Index')+6:fname.find('Index')+9])
    mag = float(fname[fname.find('Mag')+4:fname.find('Mag')+9])
    density = float(fname[fname.find('Den_')+4:fname.find('_Egal')])
    key = str(A)+'_'+str(spec)+'_'+str(density)+str(mag)
    print 'Specs ', A, Z, spec, mag, density
    arr = np.genfromtxt(fname, delimiter='|')
    if np.any(np.isnan(arr[0])):
        print 'Finding Nans'
        return np.nan
    fname = 'Outputs/'+str(A)+'_'+str(Z)+'/'+fname
    egalmap = arr[0]
    arr[1][np.isnan(arr[1])] = np.nanmean(arr[1])
    zmap = inv_comoving_distance(arr[1])
    comp5map = np.round(comp5dict[key](zmap))
    comp1map = np.round(comp1dict[key](zmap))
    comp9map = np.round(comp9dict[key](zmap))
    compsdmap = (comp9map-comp1map)/4.
    if np.any(np.isnan(compsdmap)):
        print 'Nans in compsdmap'
    if len(compsdmap[compsdmap==0.]):
        print 'Zeros in compsdmap'
    comp5map = (comp5map+Z)/2.
    
    
    
    print 'Min, Med, Max of composition', np.min(comp1map), np.median(comp5map), np.max(comp9map)
    mincomp, medcom, maxcomp = np.min(comp1map), np.median(comp5map), np.max(comp9map)
    if mincomp>Z:
        mincomp = Z-10
    reqcomps=range(int(mincomp), int(Z+1))
    en5map = en5dict[key](zmap)
    en1map = en1dict[key](zmap)
    en9map = en9dict[key](zmap)
    print 'Min, Med, Max of Energy', np.min(en1map), np.median(en5map), np.max(en9map)
    
    
    if not len(evdict.keys()):
        print 'Looking for all Galactic backtrack files in the z range:', reqcomps[0], reqcomps[-1]
        for comp in reqcomps:
            print 'Now looking for composition, ', comp
            print reqcomps
            if not comp in nucmassdict.keys():
                continue
            fwc = 'Outputs/'+str(nucmassdict[int(comp)])+'_'+str(int(comp))+'/MinEn53.0_N10_SpecI4.3_Seq7*GalBack_Manual.txt'    
            flist = sorted(glob(fwc))
            print 'found ', len(flist), 'files for composition ', comp
            if not len(flist):
                continue
            try:
                evarr = np.genfromtxt(flist[0], delimiter="|").transpose()
                if not len(evarr):
                    continue
                for fnamen in flist[1:]:
                    print fnamen
                    now = np.genfromtxt(fnamen, delimiter="|").transpose()
                    print len(now[0]), 'events found'
                    evarr = np.append(evarr, now , axis=1)
                evdict[comp] = evarr
                evweightdict[comp] = 1./float(len(flist))
            except:
                print 'Bah'

        
    posweightedmap = {}
    print 'Composition weighting now'
    mainmap = np.zeros(hp.nside2npix(32))
    pltcount=0
    plt.figure(pltcount)
    hp.mollview(egalmap, title='Map outside Galaxy')
    plt.savefig(fname.replace('.txt', '.png'))
    for comp in sorted(evdict.keys()):
            evarr = evdict[comp]
            coord = SkyCoord(l = evarr[1]*u.radian, b = evarr[0]*u.radian, frame='galactic')
            ra = coord.icrs.ra.value
            dec = coord.icrs.dec.value

            defcoord = SkyCoord(l = evarr[4]*u.radian, b = evarr[3]*u.radian, frame='galactic')
            defra = defcoord.icrs.ra.value
            defdec = defcoord.icrs.dec.value

            defmap = scattomap(defdec, defra, nside=32)
            #plt.figure(pltcount)
            #hp.mollview(defmap, title='Backpropagated Map at Galactic border, isotropic at Earth')
            #plt.show()
            #plt.savefig(fname.replace('.txt', str(comp)+'JusDefs.png'))
            #pltcount=pltcount+1
            defpixs = hp.ang2pix(32, np.deg2rad(90.-defdec), np.deg2rad(defra))
            compweight = 1./(compsdmap*np.sqrt(2.*np.pi))*np.exp(-0.5*(comp5map - comp)*(comp5map - comp)/(compsdmap*compsdmap))
            if np.any(np.isnan(compweight)):
                print 'Nans in compweight'
                plt.figure(pltcount)
                hp.mollview(compweight, title=str(comp))
                plt.show()
                pltcount+=1
                
                compweight[np.isnan(compweight)]=np.nanmean(compweight)
            if A==1:
                compweight = np.ones(hp.nside2npix(32))
            weights = egalmap[defpixs]#Add Spectral weight map here
            posweightedmap[comp] = weightedscattomap(dec, ra, weights, 32)*evweightdict[comp]*compweight
            if np.any(np.isnan(posweightedmap[comp])):
                print comp, evweightdict[comp], compweight
                
            mainmap = mainmap+posweightedmap[comp]
    
    mainmap = mainmap/np.sum(mainmap)
    
    print mainmap, len(mainmap)
    
    plt.figure(pltcount)
    hp.mollview(mainmap, title='Earth Map')
    np.savetxt(fname.replace('EGal_Injection_', 'AtEarthMap2_'), mainmap)
    plt.savefig(fname.replace('.txt', str(comp)+'EarthMap.png'))
    #plt.show()
    del plt
    return mainmap
    
if len(allflist):
    fout = open('Outputs/GalSummary5_'+str(A)+'_'+str(Z)+'_SpecIndex_'+str(spectralindex)+'_Den_'+str(Density)+'_EgalMag_'+str(magfield)+'.txt', "w")

for toproc in allflist:
    mapinstance = ProcessMap(toproc)
    print np.sum(mapinstance), mapinstance
    hp.mollview(mapinstance)
    plt.show()
    try:
        if not np.any(np.isnan(mapinstance)):
            llh = logLikelihood(mapinstance)
            print toproc, llh
        else:
            llh = np.nan
        fout.write(toproc+ '|' + str(llh)+ ' \n')
    except:
        print 'Something Fucked up'
    
fout.close()
    


"""
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

"""
