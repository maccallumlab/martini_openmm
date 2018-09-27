An implementation of the Martini CG force field in OpenMM.

This is alpha quality software. Use at your own risk. Testing
on your specific system is recommend. See the `tests` directory
for example scripts to compare gromacs and openmm energies and
forces.

Please report any defects as issues on
[the github repo](https://github.com/maccallumlab/martini_openmm).

The system has been tested on the sytems in the `tests` directory
and works for:
- Martini v2
- Simple and complex mixtures of lipids
- Simple simulations of soluble proteins

It has not yet been tested for:
- Elastic network models of proteins
- ElNeDyn models of proteins

There has not yet been any effort to support:
- Martini v3
- Polarizable Martini