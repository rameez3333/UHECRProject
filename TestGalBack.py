from crpropa import *
import os
from optparse import OptionParser
import healpy as hp
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord


usage = 'usage: %prog [options]'
parser = OptionParser(usage)

parser.add_option("-a", "--massnumber", action="store", type="int", default=1, dest="MASSNUMBER", help="The Atomic Mass Number (A) of the Nucleus")
parser.add_option("-z", "--atomicnumber", action="store", type="int", default=1, dest="ATOMICNUMBER", help="The Atomic Number (Z) of the Nucleus")
parser.add_option("-m", "--minenerg", action="store", type="float", default=53.0, dest="MINENERG", help="The minimum energy of the CR beyond which it is not tracked (in EeV)")
parser.add_option("-n", "--howmany", action="store", type="int", default=1, dest="HOWMANY", help="How many CRs to propagate")
parser.add_option("-f", "--firstseq", action="store", type="int", default=0, dest="FIRSTSEQ", help="The sequence number")

(options, args) = parser.parse_args()


massno = options.MASSNUMBER
atomicno = options.ATOMICNUMBER
minenerg = options.MINENERG
howmany = options.HOWMANY
sequence = options.FIRSTSEQ

def generatepowerlaw(index, rmin, rmax, size):
    arr = np.random.uniform(np.power(rmin, -1.*index), np.power(rmax, -1.*index), size)
    return np.power(arr, -1./index)


# magnetic field setup
B = JF12Field()
seed = 691342
B.randomStriated(seed)
B.randomTurbulent(seed)

# simulation setup
sim = ModuleList()
sim.add(PropagationCK(B, 1e-4, 0.1 * parsec, 100 * parsec))
obs = Observer()
obs.add(ObserverLargeSphere(Vector3d(0), 20 * kpc))
# obs.onDetection(TextOutput('galactic_backtracking.txt', Output.Event3D))
sim.add(obs)
print sim

nside=128

output = TextOutput('Outputs/'+str(massno)+'_'+str(atomicno)+'/MinEn'+str(minenerg)+'_N'+str(howmany)+'_Seq'+str(sequence)+'GalBack.txt', Output.Event3D)
sim.add(output)

detectedspectralindex = 3
energies = generatepowerlaw(detectedspectralindex, minenerg, 200., 100000) 

def propagate(A, Z, energy, lat, lon, seq=0):
    print lon, lat, energy
    pid = - nucleusId(A,Z)  # (anti-)proton
    position = Vector3d(-8.5, 0, 0) * kpc
    direction = Vector3d()
    direction.setRThetaPhi(1, lat, lon)
    p = ParticleState(pid, energy * EeV, position, direction)
    c = Candidate(p)
    sim.run(c)
    print c
    d1 = c.current.getDirection()  # direction at galactic border
    print 'galactic deflection %.2f radian' % direction.getAngleTo(d1)
    
for i in range(0, hp.nside2npix(nside)):
    dec, ra = np.deg2rad((90. - np.rad2deg(hp.pix2ang(nside, i)[0]))), hp.pix2ang(nside, i)[1]
    lat = np.deg2rad(SkyCoord(ra = ra*u.radian, dec=dec*u.radian).galactic.b.value)
    lon = np.deg2rad(SkyCoord(ra = ra*u.radian, dec=dec*u.radian).galactic.l.value)
    for i in range(howmany):
        energy = np.random.choice(energies)
        propagate(massno, atomicno, energy, lat, lon)
