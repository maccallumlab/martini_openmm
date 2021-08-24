#!/usr/bin/env python

import numpy as np
import sys


# tolerances for energy and force comparisons
energy_atol = 1e-4
energy_rtol = 1e-5
f_atol = 1e-4
f_rtol = 1e-5


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
    print(f"    Values differ at {n_diff} of {gmx_forces.shape[0]} positions.")
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
    print("Test failed!")
    sys.exit(1)
else:
    print("Test succeeded!")
