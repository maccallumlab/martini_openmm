from simtk import unit as u
from simtk import openmm as mm
from simtk.openmm import app
import martini_openmm as martini
import numpy as np

# do all tests with reference platform to avoid messing
# with cuda / drivers / etc
platform = mm.Platform.getPlatformByName("Reference")

conf = app.GromacsGroFile("minimized.gro")
box_vectors = conf.getPeriodicBoxVectors()
top = martini.GromacsMartiniV2TopFile("dppc.top", periodicBoxVectors=box_vectors)

system = top.createSystem(nonbondedCutoff=1.1 * u.nanometer)
integrator = mm.LangevinIntegrator(
    300 * u.kelvin, 1.0 / u.picosecond, 2 * u.femtosecond
)
simulation = mm.app.Simulation(top.topology, system, integrator, platform)

simulation.context.setPositions(conf.getPositions())
state = simulation.context.getState(getEnergy=True, getForces=True)

energy = state.getPotentialEnergy().value_in_unit(u.kilojoule_per_mole)
forces = np.array(state.getForces().value_in_unit(u.kilojoule / u.nanometer / u.mole))

with open('energy.txt', 'w') as efile:
    print(energy, file=efile)
np.savetxt('forces.txt', forces)