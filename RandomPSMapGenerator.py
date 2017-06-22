import numpy as np
import matplotlib.pyplot as plt
import healpy as hp
from glob import glob
from astropy import units as u
from astropy.coordinates import SkyCoord
from optparse import OptionParser
from scipy.optimize import fminbound
from scipy.optimize import minimize	
from scipy.interpolate import interp1d
from scipy.interpolate import InterpolatedUnivariateSpline
from astropy.cosmology import Planck15 as cosmo

usage = 'usage: %prog [options]'
parser = OptionParser(usage)

parser.add_option("-a", "--massnumber", action="store", type="int", default=1, dest="MASSNUMBER", help="The Atomic Mass Number (A) of the Nucleus")
parser.add_option("-z", "--atomicnumber", action="store", type="int", default=1, dest="ATOMICNUMBER", help="The Atomic Number (Z) of the Nucleus")
parser.add_option("-p", "--psnumber", action="store", type="int", default=30, dest="PSNUMBER", help="The Number of Point Sources")
parser.add_option("-u", "--unumber", action="store", type="int", default=300, dest="UNUMBER", help="The Number of UHECR events")
parser.add_option("-s", "--seqnumber", action="store", type="int", default=1, dest="SEQ", help="Sequence for output file")
parser.add_option("-e", "--ethresh", action="store", type="float", default=57.0, dest="ETHRESH", help="Sequence for output file")
parser.add_option("-w", "--weightfordist", action = "store_true", default=False, dest="WEIGHT", help = "Weight the point sources according to an evolution model and propagation effects?")


(options, args) = parser.parse_args()

massno = options.MASSNUMBER
atomicno = options.ATOMICNUMBER
psno = options.PSNUMBER
uhecrno = options.UNUMBER
seqno = options.SEQ
ethresh = options.ETHRESH
weightflag = options.WEIGHT

fname = str(massno)+'_'+str(atomicno)+'_p'+str(psno)+'_u'+str(uhecrno)+'_E'+str(ethresh)+'_'+str(seqno)

def LH(x) :
	Bx = x[0]
	By = x[1]
	Bz = x[2]
	fraction = x[3]
	thetascale = np.exp(x[4])
	PSscale = np.exp(scalefunclog(np.log(0.5)))
	
	B0 = np.sqrt(Bx*Bx + By*By + Bz*Bz)
	#B0 = 0.0
	
	if B0 > 0.0 :
	
		Bx = Bx/B0
		By = By/B0
		Bz = Bz/B0
	
		proj = Bx*vx + By*vy + Bz*vz
		
		vxtemp1 = proj*Bx
		vytemp1 = proj*By
		vztemp1 = proj*Bz
	
		vxtemp2 = vx - vxtemp1
		vytemp2 = vy - vytemp1
		vztemp2 = vz - vztemp1
	
		vxtemp3 = By*vz - Bz*vy
		vytemp3 = Bz*vx - Bx*vz
		vztemp3 = Bx*vy - By*vx
	
	total = 0.0
	
	for j in range(0,len(VXtemp)) :
			
		theta50 = thetascale*(100.0/ENtemp[j]) # rescaled with rigidity
		
		if theta50 >= 90.0 : # maximum theta50 scale
			scale = 0.0
		elif theta50 <= 0.5 : # minimum theta50 scale for interpolation
			scale = np.exp(scalefunclog(np.log(0.5)))
		elif theta50 > 60.0 : # lin-lin interpolation better above 60deg
			scale = scalefunclin(theta50)
		else :  # else log-log interpolation
			scale = np.exp(scalefunclog(np.log(theta50)))
			
		# account for angular resolution:	
		scale = scale/(1.0 + scale/PSscale)
		
		wx = VXtemp[j]
		wy = VYtemp[j]
		wz = VZtemp[j]
		
		if B0 > 0.0 :
			phase = B0*(100.0/ENtemp[j])
		
			temp1 =  vxtemp1*wx + vytemp1*wy + vztemp1*wz
			temp2 = (vxtemp2*wx + vytemp2*wy + vztemp2*wz)*np.cos(phase)
			temp3 = (vxtemp3*wx + vytemp3*wy + vztemp3*wz)*np.sin(phase)
			
			temp = temp1 + temp2 + temp3
		else :
			temp = vx*wx + vy*wy + vz*wz
			
		#calculate S/B for events:
		if scale > 0.01 :
			SoverB = 2.0*scale/(1.0-np.exp(-2.*scale))*np.exp((temp-1.0)*scale)
		else : # first order term for scale<<1
			SoverB = np.exp((temp-1.0)*scale)
			
		prob = fraction*SoverB + (1.0-fraction)
		
		if fraction < 1.0 and prob > 0.0 :
			total += np.log(prob)
		
	# correct for the limited event number after angle cut:	
	if fraction < 1.0 :	
		total += (len(VX)*1.0-len(VXtemp)*1.0)*np.log(1.0-fraction)
	
	# add Gaussian priors to regular and random fields	
	priors = (B0*180./np.pi/10.0)**2 + (thetascale/10.0)**2
	
	# return LLH:
	return -2.0*total + priors

def GetRandomPosition(size=1):
    ra = np.random.uniform(0, 360, size=size)
    dec = np.rad2deg(np.arcsin(np.random.uniform(0, 1., size=size))*np.random.choice([1.,-1.], size=size))
    if size==1:
        ra = ra[0]
        dec = dec[0]
    
    return dec,ra

def scattomap(dec,ra, nside=32):
    hmap = np.bincount(hp.ang2pix(nside, np.deg2rad(90.-dec), np.deg2rad(ra)), minlength=hp.nside2npix(nside))
    return hmap


def weightedscattomap(dec,ra, weights, nside=32):
    hmap = np.bincount(hp.ang2pix(nside, np.deg2rad(90.-dec), np.deg2rad(ra)), weights=weights, minlength=hp.nside2npix(nside))
    return hmap    

def plot_mwd(RA,Dec,org=0,title='Mollweide projection', projection='mollweide'):
    ''' RA, Dec are arrays of the same length.                                                                                                                                                                                                          
    RA takes values in [0,360), Dec in [-90,90],                                                                                                                                                                                                        
    which represent angles in degrees.                                                                                                                                                                                                                  
    org is the origin of the plot, 0 or a multiple of 30 degrees in [0,360).                                                                                                                                                                            
    title is the title of the figure.                                                                                                                                                                                                                   
    projection is the kind of projection: 'mollweide', 'aitoff', 'hammer', 'lambert'                                                                                                                                                                    
    '''
    x = np.remainder(RA+360-org,360) # shift RA values                                                                                                                                                                                                  
    ind = x>180
    x[ind] -=360    # scale conversion to [-180, 180]                                                                                                                                                                                                   
    x=-x    # reverse the scale: East to the left                                                                                                                                                                                                       
    tick_labels = np.array([150, 120, 90, 60, 30, 0, 330, 300, 270, 240, 210])
    tick_labels = np.remainder(tick_labels+360+org,360)
    fig = plt.figure(figsize=(10, 5))
    ax = fig.add_subplot(111, projection=projection, axisbg ='LightCyan')
    ax.scatter(np.radians(x),np.radians(Dec), alpha=0.8, s=0.5)  # convert degrees to radians                                                                                                                                                           
    ax.set_xticklabels(tick_labels)     # we add the scale on the x axis                                                                                                                                                                                
    ax.set_title(title)
    ax.title.set_fontsize(15)
    ax.set_xlabel("RA")
    ax.xaxis.label.set_fontsize(12)
    ax.set_ylabel("Dec")
    ax.yaxis.label.set_fontsize(12)
    ax.grid(True, which='both')

psdecs, psras = GetRandomPosition(psno)

if weightflag:
    print 'Generating according to evolution model'
    #Generate redshifts according toavolution
    M=3
    zmax = 0.13
    Evoprobx = np.linspace(0, zmax, 1500)
    Evoproby = np.power(1.+Evoprobx, M)*np.power(Evoprobx, 2.)
    backpolate = InterpolatedUnivariateSpline(Evoproby, Evoprobx)
    rs = np.random.uniform(0, np.power((1.+zmax), M)*np.power(zmax,2), size=psno)
    z = backpolate(rs)
    dists=cosmo.comoving_distance(z).value
    proparr = np.genfromtxt('1DPropagationProbability_1_1_4.3_2.txt', delimiter='|')
    propprob = InterpolatedUnivariateSpline(proparr[1], proparr[0], ext=3)
    geoweight = 1./(4.*np.pi*dists*dists)
    propweight = propprob(dists)
    totweight = geoweight*propweight
    print 'weights', totweight
    print z
    print dists
    print geoweight
    print propweight
    psmap = weightedscattomap(psdecs, psras, totweight)
else:
    psmap = scattomap(psdecs, psras)

psmap = hp.smoothing(psmap, sigma=np.deg2rad(2.0), verbose=False)
psmap = psmap+-1.*np.min(psmap)
psmap = psmap/np.sum(psmap)

flist = sorted(glob('Outputs/1_1/MinEn53.0_N10_SpecI4.3_Seq70?GalBack_Manual.txt'))

print len(flist), 'files found'

arr = np.genfromtxt(flist[0], delimiter="|").transpose()
print len(arr[0]), 'events found'


for finame in flist[1:]:
    print finame
    now = np.genfromtxt(finame, delimiter="|").transpose()
    print len(now[0]), 'events found'
    arr = np.append(arr, now , axis=1)
    
coord = SkyCoord(l = arr[1]*u.radian, b = arr[0]*u.radian, frame='galactic')
ra = coord.icrs.ra.value
dec = coord.icrs.dec.value

defcoord = SkyCoord(l = arr[4]*u.radian, b = arr[3]*u.radian, frame='galactic')

defra = defcoord.icrs.ra.value
defdec = defcoord.icrs.dec.value    


defpixs = hp.ang2pix(32, np.deg2rad(90.-defdec), np.deg2rad(defra))

weights = psmap[defpixs]


weights = weights/float(np.sum(weights))
pick = np.random.choice(len(arr[0]), uhecrno, replace=False, p=weights)

pickedarr = arr.transpose()[pick]

pickedarr = pickedarr.transpose()

print 'GeneratedUHECRs:', pickedarr
print 'Generated Point Sources :', psdecs, psras

pickedcoord = SkyCoord(l = pickedarr[1]*u.radian, b = pickedarr[0]*u.radian, frame='galactic')
pickedra = pickedcoord.icrs.ra.value
pickeddec = pickedcoord.icrs.dec.value   

plot_mwd(psras, psdecs, title='Point Sources')

#plt.show()
plt.savefig('PS_'+fname+'.png')
if weightflag:
    np.savetxt('PS_'+fname+'.txt', np.vstack([psras, psdecs, z]).transpose(), delimiter='|')
else:
    np.savetxt('PS_'+fname+'.txt', np.vstack([psras, psdecs]).transpose(), delimiter='|')
del plt

import matplotlib.pyplot as plt

plot_mwd(pickedra, pickeddec, title = 'UHECRs')
#plt.show()
plt.savefig('UHECR_'+fname+'.png')
np.savetxt('UHECR_'+fname+'.txt', pickedarr.transpose(), delimiter='|')


lines = [line.strip() for line in open("theta_to_scale.dat")]
array = [line.split("\t") for line in lines]
array2 = [[float(element) for element in line] for line in array]
A = np.array(array2)
B0 = A.transpose()

scalefunclog = interp1d(np.log(B0[0][:len(B0[0])-1]), np.log(B0[1][:len(B0[0])-1]))
scalefunclin = interp1d(B0[0], B0[1])

nside = 32

npix = hp.nside2npix(nside)

map = np.zeros(npix, dtype=np.double)
mapx = np.zeros(npix, dtype=np.double)
mapy = np.zeros(npix, dtype=np.double)
mapz = np.zeros(npix, dtype=np.double)
mapns = np.zeros(npix, dtype=np.double)
maptheta = np.zeros(npix, dtype=np.double)

Ecut = ethresh
ANGLECUT = 50.0 
pickedarr = pickedarr.transpose()[pickedarr[2]>Ecut].transpose()

pickedcoord = SkyCoord(l = pickedarr[1]*u.radian, b = pickedarr[0]*u.radian, frame='galactic')
RA = pickedcoord.icrs.ra.value
DEC = pickedcoord.icrs.dec.value   
EN = pickedarr[2]
temp = hp.ang2vec(np.deg2rad(90.-DEC), np.deg2rad(RA)).transpose()
VX, VY, VZ = temp[0], temp[1], temp[2]

COLOR='green'
SYMBOL='x'

for i in range(0,npix) :
    
    print(i)
    
    vx,vy,vz = hp.pix2vec(nside,i)
    
    RAtemp = []
    DECtemp = []

    VXtemp = []
    VYtemp = []
    VZtemp = []

    ENtemp = []
    
    for j in range(0,len(VX)) :
    
        x = vx*VX[j] + vy*VY[j] + vz*VZ[j]
        
        angdis = np.arccos(x)/np.pi*180.0
        
        if angdis < ANGLECUT : # introduce effective data with angle cut 
            VXtemp.append(VX[j])
            VYtemp.append(VY[j])
            VZtemp.append(VZ[j])
            ENtemp.append(EN[j])
    
    start = (0.0,0.0,0.0,0.5,np.log(10.0))
    boundary = ((-1.0,1.0),(-1.0,1.0),(-1.0,1.0),(0.0,1.0),(np.log(0.1),np.log(90.0)))
    
    if len(VXtemp) > 0 :	
        res = minimize(LH, start,bounds=boundary)
        
        map[i] = -LH(res.x)
        mapns[i] = res.x[3]*len(VX) # signal events
        
        mapx[i] = res.x[0]/np.pi*180.0 # units of degree
        mapy[i] = res.x[1]/np.pi*180.0 # units of degree
        mapz[i] = res.x[2]/np.pi*180.0 # units of degree
        maptheta[i] = np.exp(res.x[4]) # units of degree

        
#hp.write_map("map_LH_"+fname+"_nside32.fits", map) # likelihood map
#hp.write_map("map_ns_"+fname+"_nside32.fits",mapns)	# signal events
#hp.write_map("map_theta_"+fname+"_nside32.fits",maptheta)	# von Mises theta50 at 100EV (in degree)
#hp.write_map("map_Bx_"+fname+"_nside32.fits",mapx) # rotation at 100EV (in degree)
#hp.write_map("map_By_"+fname+"_nside32.fits",mapy) # rotation at 100EV (in degree)
#hp.write_map("map_Bz_"+fname+"_nside32.fits",mapz) # rotation at 100EV (in degree)

np.savetxt('maps_LH'+fname+'nside32.txt', np.vstack([map, mapns, maptheta, mapx, mapy, mapz]), delimiter='|')
