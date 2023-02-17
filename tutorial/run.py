from openmm.unit import *
from openmm import *
from openmm.app import *
import martini_openmm as martini
from mdtraj.reporters import XTCReporter
from sys import stdout

def run(epsilon_r):

	platform = Platform.getPlatformByName("CUDA")
	properties = {'Precision': 'double'}

	conf = GromacsGroFile("system.gro")
	box_vectors = conf.getPeriodicBoxVectors()

	# get any defines
	defines = {}
	try:
		with open("defines.txt") as def_file:
			for line in def_file:
				line = line.strip()
				defines[line] = True
	except FileNotFoundError:
		pass

	top = martini.MartiniTopFile(
		"system.top",
		periodicBoxVectors=box_vectors,
		defines=defines,
		epsilon_r=epsilon_r,
	)
	system = top.create_system(nonbonded_cutoff=1.1 * nanometer)

	integrator = LangevinIntegrator(310 * kelvin,
									10.0 / picosecond,
									20 * femtosecond)
	integrator.setRandomNumberSeed(0)

	simulation = Simulation(top.topology, system, integrator,
							platform, properties)

	simulation.context.setPositions(conf.getPositions())

################################################################################
		### Minimization ###

	simulation.reporters.append(PDBReporter('mini.pdb', 1000))
	simulation.reporters.append(StateDataReporter(stdout, 5000,
													step=True,
													potentialEnergy=True,
													temperature=True,
													volume=True)
								)
	print("Minimizing energy...")
	simulation.minimizeEnergy(maxIterations=5000,tolerance=1.0)

	energies = simulation.context.getState(getEnergy=True).getPotentialEnergy()
	print("System minimized at", energies, "\n")

################################################################################
		### NVT equilibration ###

	simulation.context.setVelocitiesToTemperature(310 * kelvin)
	print('Running NVT equilibration...')
	simulation.step(50000) #1ns

################################################################################
		### NPT equilibration ###
	
	system.addForce(MonteCarloBarostat(1 * bar, 310 * kelvin))
	# to update the simulation object to take in account the new system
	simulation.context.reinitialize(True)
	print('Running NPT equilibration...')
	simulation.step(50000) #1ns

	# save the equilibration results to file
	simulation.saveState('equi.state')
	simulation.saveCheckpoint('equi.chk')

################################################################################
		### Production run ###

	# Set up the reporters to report energies every 1000 steps.
	simulation.reporters.append(StateDataReporter("prod.log", 1000,
													step=True,
													potentialEnergy=True,
													totalEnergy=True,
													density=True,
													temperature=True,
													volume=True)
								)
	# save the trajectory in XTC format
	xtc_reporter = XTCReporter('prod.xtc', 1000)
	simulation.reporters.append(xtc_reporter)

	# run simulation
	print("Running simulation...")
	simulation.step(50000) #1ns

run(15)
