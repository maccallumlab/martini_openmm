An implementation of the Martini CG force field in OpenMM.

This is alpha quality software. Use at your own risk. Testing
on your specific system is recommend. See the `tests` directory
for example scripts to compare gromacs and openmm energies and
forces.

A suite of tests is run automatically. The current status of these tests is:
[![martini_openmm](https://github.com/maccallumlab/martini_openmm/actions/workflows/CI.yml/badge.svg)](https://github.com/maccallumlab/martini_openmm/actions)

Please report any defects as issues on
[the github repo](https://github.com/maccallumlab/martini_openmm).

The system has been tested on the sytems in the `tests` directory
and works for:
- Martini v2
- Martini v3
- Polarizable Martini
- Simple and complex mixtures of lipids
- Simple simulations of soluble and membrane proteins
To run the tests, use `pytest` in the root directory after installation.

Limitations:
- Math in .itp files:
    - Gromacs allows for mathmatical expressions to be used in .itp files
        - `    1     6      1   0.98112 RUBBER_FC*1.000000`
    - This is not allowed in OpenMM and must be edited
        - `    1     6      1   0.98112 RUBBER_FC`
- The only supported electrostatic option is `coulombtype = reaction-field`
- The only supported option for Van der Waals is `vdw_type = cutoff`
- The cutoff must be the same for `rcoulomb` and `rvdw`
