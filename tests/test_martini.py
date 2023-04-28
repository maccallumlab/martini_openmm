from turtle import pos
import unittest
import os
import subprocess
from simtk import unit as u
from simtk import openmm as mm
from simtk.openmm import app
import martini_openmm as martini
import numpy as np


# tolerances for energy and force comparisons
energy_atol = 1e-4
energy_rtol = 1e-5
f_atol = 1e-4
f_rtol = 1e-5


class Base:
    epsilon_r = 15.0

    def setUp(self):
        self.base_dir = os.getcwd()
        os.chdir(self.test_dir)

    def tearDown(self):
        os.chdir(self.base_dir)

    def test(self):
        # Run gromacs
        prefix = "../" * self.depth
        result = subprocess.run(
            f"{prefix}Scripts/run_gmx.sh", shell=True, capture_output=True, text=True
        )
        _handle_gmx_result(result)

        # Run openmm
        _gen_openmm(self.epsilon_r)

        # Do comparison
        _compare()

        # Clean up
        _cleanup()


class TestAdenine(Base, unittest.TestCase):
    test_dir = "tests/aden_m3"
    depth = 1


class TestThymine(Base, unittest.TestCase):
    test_dir = "tests/thym_m3"
    depth = 1


class TestThymine2(Base, unittest.TestCase):
    test_dir = "tests/thy2_m3"
    depth = 1


class TestPolWater(Base, unittest.TestCase):
    test_dir = "tests/pol_water"
    depth = 1
    epsilon_r = 2.5


class TestMOL1(Base, unittest.TestCase):
    test_dir = "tests/MOL1"
    depth = 1


class TestPEG(Base, unittest.TestCase):
    test_dir = "tests/PEG"
    depth = 1


class TestMartini2SimpleLipid(Base, unittest.TestCase):
    test_dir = "tests/simple_lipid"
    depth = 1


class TestMartini2Protein(Base, unittest.TestCase):
    test_dir = "tests/protein"
    depth = 1


class TestMartini2ComplexLipid(Base, unittest.TestCase):
    test_dir = "tests/complex_lipid"
    depth = 1


class TestMartini2Cholesterol(Base, unittest.TestCase):
    test_dir = "tests/m2_chol"
    depth = 1


class TestMartini2Cholesterol(Base, unittest.TestCase):
    test_dir = "tests/m3_chol"
    depth = 1


class TestMartini2MembraneProtein(Base, unittest.TestCase):
    test_dir = "tests/membrane_protein"
    depth = 1


class TestMartini3AqpEN(Base, unittest.TestCase):
    test_dir = "tests/1j4n_en_m3"
    depth = 1


class TestMartini3UbqEN(Base, unittest.TestCase):
    test_dir = "tests/1ubq_en_m3"
    depth = 1


class TestMartini3UbqGo(Base, unittest.TestCase):
    test_dir = "tests/1ubq_go_m3"
    depth = 1


class TestMartini3BptiEN(Base, unittest.TestCase):
    test_dir = "tests/1k6u_en_m3"
    depth = 1


class TestMartini3BptiGo(Base, unittest.TestCase):
    test_dir = "tests/1k6u_go_m3"
    depth = 1


class TestMartini3POPC(Base, unittest.TestCase):
    test_dir = "tests/popc_m3"
    depth = 1


class TestMartini3Dodecane(Base, unittest.TestCase):
    test_dir = "tests/others_m3/1DOD"
    depth = 2


class TestMartini3NACL(Base, unittest.TestCase):
    test_dir = "tests/others_m3/1NACL"
    depth = 2


class TestMartini3SinglePOPC(Base, unittest.TestCase):
    test_dir = "tests/others_m3/1POPC"
    depth = 2


class TestMartini3Water(Base, unittest.TestCase):
    test_dir = "tests/others_m3/2W"
    depth = 2


class TestMartini3SmallMono3HT(Base, unittest.TestCase):
    test_dir = "tests/small_mols_m3/mono/3HT"
    depth = 3


class TestMartini3SmallMono4NIAN(Base, unittest.TestCase):
    test_dir = "tests/small_mols_m3/mono/4NIAN"
    depth = 3


class TestMartini3SmallMonoCLPR(Base, unittest.TestCase):
    test_dir = "tests/small_mols_m3/mono/CLPR"
    depth = 3


class TestMartini3SmallPoly2T(Base, unittest.TestCase):
    test_dir = "tests/small_mols_m3/poly/2T"
    depth = 3


class TestMartini3SmallPolyANTH(Base, unittest.TestCase):
    test_dir = "tests/small_mols_m3/poly/ANTH"
    depth = 3


class TestMartini3SmallPolyCAFF(Base, unittest.TestCase):
    test_dir = "tests/small_mols_m3/poly/CAFF"
    depth = 3


class TestMartini3SmallPolyCNAP(Base, unittest.TestCase):
    test_dir = "tests/small_mols_m3/poly/CNAP"
    depth = 3


class TestMartini3SmallPolyNDMBI(Base, unittest.TestCase):
    test_dir = "tests/small_mols_m3/poly/NDMBI"
    depth = 3


class TestMartini3SmallPolyTDMBI(Base, unittest.TestCase):
    test_dir = "tests/small_mols_m3/poly/TDMBI"
    depth = 3

class TestMartini2ElasticNetworkType6(Base, unittest.TestCase):
    test_dir = "tests/1ubq_en_m2_btype6"
    depth = 1


def _handle_gmx_result(result):
    if result.returncode != 0:
        print("Gromacs failed to run")
        print()
        print()
        print("Standard Output:")
        print()
        print(result.stdout)
        print()
        print("Standard Error")
        print()
        print(result.stderr)
        print()
        raise RuntimeError(
            "Gromacs failed to run. See standard output and error above."
        )


def _gen_openmm(epsilon_r):
    os.chdir("openmm")

    # do all tests with reference platform to avoid messing
    # with cuda / drivers / etc
    # this also uses full double precision
    platform = mm.Platform.getPlatformByName("Reference")

    conf = app.GromacsGroFile("minimized.gro")
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
    system = top.create_system(nonbonded_cutoff=1.1 * u.nanometer)
    integrator = mm.LangevinIntegrator(
        300 * u.kelvin, 1.0 / u.picosecond, 2 * u.femtosecond
    )

    simulation = mm.app.Simulation(top.topology, system, integrator, platform)

    # Check to see that constraints and vsites are applied correctly
    simulation.context.setPositions(conf.getPositions())
    simulation.context.applyConstraints(tol=1e-5)
    state = simulation.context.getState(getPositions=True)
    positions = state.getPositions(asNumpy=True).value_in_unit(u.nanometer)
    ref = np.array(conf.getPositions().value_in_unit(u.nanometer))
    delta = np.linalg.norm(positions - ref, axis=1)
    if np.any(delta > 2.5e-3):
        index = np.argmax(delta)
        start = max(0, index - 5)
        end = min(positions.shape[0], index + 8)
        print(positions[start:end, :])
        raise AssertionError(
            f"""
            Applying constraints or vsites resulted in a position error
            of > 2.5e-3 nm. Largest deviation is on atom {index}
            (zero-based) with deviation {delta[index]}.
            """
        )

    # Compute energies and forces
    simulation.context.setPositions(conf.getPositions())
    state = simulation.context.getState(getEnergy=True, getForces=True)
    energy = state.getPotentialEnergy().value_in_unit(u.kilojoule_per_mole)
    forces = np.array(
        state.getForces().value_in_unit(u.kilojoule / u.nanometer / u.mole)
    )

    # Openmm does not integrate the positions of the virtual sites. It updates them
    # after intetrating the positions of other atoms. As such, it does not set the
    # forces on virtual sites to zero. This has no effect on the trajectory, but
    # GROMACS reports these forces as zero. To check for consistency, we therefor
    # need to set these forces to zero.
    for site in top._all_vsites:
        forces[site, :] = 0.0

    with open("energy.txt", "w") as efile:
        print(energy, file=efile)
    np.savetxt("forces.txt", forces)

    os.chdir("..")


def _compare():
    gmx_energy = get_gmx_energy()
    omm_energy = get_omm_energy()
    gmx_forces = get_gmx_forces()
    omm_forces = get_omm_forces()

    failed = False
    relative_energy_error = abs(gmx_energy - omm_energy) / abs(gmx_energy)
    abs_energy_error = abs(gmx_energy - omm_energy)
    if relative_energy_error > 1e-5 and abs_energy_error > 1e-3:
        print("Gromacs and OpenMM energies do not match!")
        print(f"    Gromacs: {gmx_energy:16.3f}")
        print(f"     OpenMM: {omm_energy:16.3f}")
        print(f"      Delta: {(gmx_energy - omm_energy):16.3f}")
        print(f"   Relative: {relative_energy_error:16.3f}")
        print()
        failed = True

    if not np.allclose(omm_forces, gmx_forces, rtol=f_rtol, atol=f_atol):
        errors = np.logical_not(
            np.isclose(omm_forces, gmx_forces, rtol=f_rtol, atol=f_atol)
        )
        n_diff = np.sum(errors)
        bad_ind = np.where(errors)
        first_row = bad_ind[0][0]
        first_col = bad_ind[1][0]
        g_f = gmx_forces[first_row, :]
        o_f = omm_forces[first_row, :]
        d_f = np.abs(g_f - o_f)
        errors = np.abs(omm_forces - gmx_forces) - f_atol - f_rtol * np.abs(gmx_forces)
        max_ind = np.unravel_index(np.argmax(errors, axis=None), errors.shape)
        max_g = gmx_forces[max_ind[0], :]
        max_o = omm_forces[max_ind[0], :]
        max_d = np.abs(max_g - max_o)
        print("Gromacs and OpenMM forces do not match!")
        print(
            f"    Values differ at {n_diff} of {3 * gmx_forces.shape[0]} coordinates."
        )
        print()
        print(f"    First differing position: ({first_row}, {first_col})")
        print(f"    Gromacs: {g_f[0]:20f} {g_f[1]:20f} {g_f[2]:20f}")
        print(f"     OpenMM: {o_f[0]:20f} {o_f[1]:20f} {o_f[2]:20f}")
        print(f" Difference: {d_f[0]:20f} {d_f[1]:20f} {d_f[2]:20f}")
        print()
        print(f"    Largest violation is at: ({max_ind[0]}, {max_ind[1]}).")
        print(f"    Gromacs: {max_g[0]:20f} {max_g[1]:20f} {max_g[2]:20f}")
        print(f"     OpenMM: {max_o[0]:20f} {max_o[1]:20f} {max_o[2]:20f}")
        print(f" Difference: {max_d[0]:20f} {max_d[1]:20f} {max_d[2]:20f}")
        print()
        failed = True

    if failed:
        raise AssertionError("Energies or forces do not match.")


def get_gmx_energy():
    with open("gmx/energy.xvg") as f:
        line = f.readlines()[-1]
        cols = line.split()
        energy = float(cols[1])
    with open("gmx/energy.txt", "w") as outfile:
        outfile.write(f"{energy}")
    return energy


def get_omm_energy():
    with open("openmm/energy.txt") as f:
        line = f.read()
    energy = float(line)
    return energy


def get_gmx_forces():
    with open("gmx/forces.xvg") as f:
        line = f.readlines()[-1]
    cols = line.split()
    forces = np.array([float(c) for c in cols][1:])
    n = len(forces)
    forces = forces.reshape((n // 3, 3))
    np.savetxt("gmx/forces.txt", forces)
    return forces


def get_omm_forces():
    return np.loadtxt("openmm/forces.txt")


def _cleanup():
    os.remove("gmx/energy.xvg")
    os.remove("gmx/forces.xvg")
    os.remove("gmx/energy.txt")
    os.remove("gmx/forces.txt")
    os.remove("openmm/energy.txt")
    os.remove("openmm/forces.txt")
