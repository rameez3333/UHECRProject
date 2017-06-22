
from crpropa import *
import os
from optparse import OptionParser
import numpy as np



usage = 'usage: %prog [options]'
parser = OptionParser(usage)
    
parser.add_option("-a", "--massnumber", action="store", type="int", default=1, dest="MASSNUMBER", help="The Atomic Mass Number (A) of the Nucleus")
parser.add_option("-z", "--atomicnumber", action="store", type="int", default=1, dest="ATOMICNUMBER", help="The Atomic Number (Z) of the Nucleus")
parser.add_option("-m", "--minenerg", action="store", type="float", default=53.0, dest="MINENERG", help="The minimum energy of the CR beyond which it is not tracked (in EeV)")
parser.add_option("-e", "--maxenerg", action="store", type="float", default=200.0, dest="MAXENERG", help="Maximum injection CR energy")
parser.add_option("-s", "--spectralindex", action="store", type="float", default=2.0, dest="SPECTRALINDEX", help="The Spectral Index at Injection")
parser.add_option("-d", "--distance", action="store", type="float", default=100, dest="DIST", help="Starting distance of the CR (in MPc)")                  
parser.add_option("-n", "--howmany", action="store", type="int", default=1, dest="HOWMANY", help="How many CRs to propagate")
parser.add_option("-f", "--firstseq", action="store", type="int", default=0, dest="FIRSTSEQ", help="The sequence number")

(options, args) = parser.parse_args()

massno = options.MASSNUMBER
atomicno = options.ATOMICNUMBER
minenerg = options.MINENERG
maxenerg = options.MAXENERG
spectralindex = options.SPECTRALINDEX
distance = options.DIST
howmany = options.HOWMANY
sequence = options.FIRSTSEQ


detectedspectralindex = 3


sim = ModuleList()
if not os.path.isdir('Outputs/'+str(massno)+'_'+str(atomicno)):
    print 'Destination path not found. Making it now...'
    os.mkdir('Outputs/'+str(massno)+'_'+str(atomicno))
sim.add( SimplePropagation() )
sim.add( Redshift() )
# add interaction modules
sim.add( PhotoPionProduction(CMB) )
sim.add( PhotoPionProduction(IRB) )

sim.add( ElectronPairProduction(CMB) )
sim.add( ElectronPairProduction(IRB) )

sim.add( PhotoDisintegration(CMB) )
sim.add( PhotoDisintegration(IRB) )

sim.add( NuclearDecay() )

sim.add( MinimumEnergy( minenerg * EeV) )

obs = Observer()
obs.add( ObserverPoint() )  # observer at x = 0

output = TextOutput('Outputs/'+str(massno)+'_'+str(atomicno)+'/MinEn'+str(minenerg)+'MaxEn'+str(maxenerg)+'_InSpecIndex'+str(spectralindex)+'_StartDist'+str(distance)+'_N'+str(howmany)+'_Seq'+str(sequence)+'.txt', Output.Event1D)
sim.add(output)

def propagate(A, Z, startenerg, startdist, seq=0):
    cosmicray = Candidate(nucleusId(A,Z), startenerg * EeV, Vector3d(startdist * Mpc, 0, 0))
    sim.run(cosmicray)
    print cosmicray
    print 'Propagated distance', cosmicray.getTrajectoryLength() / Mpc, 'Mpc'
    
def generatepowerlaw(index, rmin, rmax, size):
    arr = np.random.uniform(np.power(rmin, -1.*index), np.power(rmax, -1.*index), size)
    return np.power(arr, -1./index)

inenergies = generatepowerlaw(spectralindex, minenerg, maxenerg, howmany)

for i in range(0, howmany):
    propagate(massno, atomicno, inenergies[i], distance, i)
