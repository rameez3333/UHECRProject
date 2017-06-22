from crpropa import *


#simulation: a sequence of simulation modules
sim = ModuleList()

# add propagator for rectalinear propagation
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

sim.add( MinimumEnergy( 90.0 * EeV) )

cosmicray = Candidate(nucleusId(1,1), 100. * EeV, Vector3d(200 * Mpc, 0, 0))

obs = Observer()
obs.add( ObserverPoint() )  # observer at x = 0

output = TextOutput('Test.txt', Output.Event1D)
#obs.onDetection( output )

sim.add(output)
print obs

sim.run(cosmicray)

print cosmicray
print 'Propagated distance', cosmicray.getTrajectoryLength() / Mpc, 'Mpc'
