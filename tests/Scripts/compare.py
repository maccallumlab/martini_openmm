#!/usr/bin/env python

import math
import numpy as np
import sys


def get_gmx_energy():
    f = open("gmx/energy.xvg")
    line = f.readlines()[-1]
    cols = line.split()
    energy = float(cols[1])
    open("gmx/energy.txt", "w").write(f"{energy}")
    return energy


def get_omm_energy():
    f = open("openmm/energy.txt")
    line = f.read()
    return float(line)


def get_gmx_forces():
    f = open("gmx/forces.xvg")
    line = f.readlines()[-1]
    cols = line.split()
    forces = np.array([float(c) for c in cols][1:])
    n = len(forces)
    forces = forces.reshape((n // 3, 3))
    np.savetxt("gmx/forces.txt", forces)
    return forces


def get_omm_forces():
    return np.loadtxt("openmm/forces.txt")


gmx_energy = get_gmx_energy()
omm_energy = get_omm_energy()
gmx_forces = get_gmx_forces()
omm_forces = get_omm_forces()
failed = False

if not math.isclose(gmx_energy, omm_energy, rel_tol=1e-6):
    print("Gromacs and OpenMM energies do not match!")
    print(f"    Gromacs: {gmx_energy}")
    print(f"     OpenMM: {omm_energy}")
    print(f"      Delta: {gmx_energy - omm_energy}")
    print()
    failed = True

if not np.allclose(gmx_forces, omm_forces, atol=2e-2):
    diffs = np.where(np.abs(gmx_forces - omm_forces) > 1e-2)
    n_diff = len(diffs[0]) // 3
    i = diffs[0][0]
    g_f = gmx_forces[i, :]
    o_f = omm_forces[i, :]
    max_diff = np.max(np.abs(gmx_forces - omm_forces))
    max_diff_ind = np.unravel_index(
        np.argmax(np.abs(gmx_forces - omm_forces)), gmx_forces.shape
    )
    ind = max_diff_ind[0]
    print("Gromacs and OpenMM forces do not match!")
    print(f"    Values differ at {n_diff} of {gmx_forces.shape[0]} positions.")
    print()
    print(f"    First differing position: {i}")
    print(f"    Gromacs: {g_f[0]:20f} {g_f[1]:20f} {g_f[2]:20f}")
    print(f"     OpenMM: {o_f[0]:20f} {o_f[1]:20f} {o_f[2]:20f}")
    print()
    print(f"    Largest difference is: {max_diff}.")
    print(f"    Largest difference is at: {ind}.")
    print(
        f"    Gromacs: {gmx_forces[ind, 0]:20f} {gmx_forces[ind, 1]:20f} {gmx_forces[ind, 2]:20f}"
    )
    print(
        f"     OpenMM: {omm_forces[ind, 0]:20f} {omm_forces[ind, 1]:20f} {omm_forces[ind, 2]:20f}"
    )
    print()
    failed = True

if failed:
    print("Test failed!")
    sys.exit(1)
else:
    print("Test succeeded!")
