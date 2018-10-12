from simtk import unit as u
from simtk import openmm as mm
from simtk.openmm import app
import martini_openmm as martini
import numpy as np


# do all tests with reference platform to avoid messing
# with cuda / drivers / etc
# this also uses full double precision
platform = mm.Platform.getPlatformByName("Reference")

conf = app.GromacsGroFile("minimized.gro")
box_vectors = conf.getPeriodicBoxVectors()

# get any defines
defines = {}
try:
    def_file = open("defines.txt")
    for line in def_file:
        line = line.strip()
        defines[line] = True
except FileNotFoundError:
    pass

top = martini.GromacsMartiniV2TopFile(
    "system.top", periodicBoxVectors=box_vectors, defines=defines
)

system = top.createSystem(nonbondedCutoff=1.1 * u.nanometer)
integrator = mm.LangevinIntegrator(
    300 * u.kelvin, 1.0 / u.picosecond, 2 * u.femtosecond
)

simulation = mm.app.Simulation(top.topology, system, integrator, platform)

simulation.context.setPositions(conf.getPositions())
state = simulation.context.getState(getEnergy=True, getForces=True)
energy = state.getPotentialEnergy().value_in_unit(u.kilojoule_per_mole)
forces = np.array(state.getForces().value_in_unit(u.kilojoule / u.nanometer / u.mole))

# Openmm does not integrate the positions of the virtual sites. It updates them
# after intetrating the positions of other atoms. As such, it does not set the
# forces on virtual sites to zero. This has no effect on the trajectory, but
# GROMACS reports these forces as zero. To check for consistency, we therefor
# need to set these forces to zero.
for atoms in top._all_vsites:
    site = atoms[0]
    forces[site, :] = 0.0

with open("energy.txt", "w") as efile:
    print(energy, file=efile)
np.savetxt("forces.txt", forces)
