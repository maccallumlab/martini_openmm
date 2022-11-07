# An implementation of the Martini coarse-grained force fields in OpenMM.

This repository describes an implementation of the Martini 2[^1][^2] and Martini 3[^3] force fields for coarse-grained simulation in the OpenMM package[^4]. It includes a Python parser to convert GROMACS Martini topology files (.top) to OpenMM topologies and systems. 

## Features
The repository includes a suite of test systems to verify the implementation of all the potential terms in use by Martini 2 and Martini 3.

The suite of tests is run automatically. The current status of these tests is:
[![martini_openmm](https://github.com/maccallumlab/martini_openmm/actions/workflows/CI.yml/badge.svg)](https://github.com/maccallumlab/martini_openmm/actions)

The test systems include:
- Simple and more complex mixtured of lipids;
- Soluble membrane proteins, with elastic networks setups for Martini 2 and 3, and including Go model for Martini 3;
Membrane proteins with elastic network setup for Martini 2 and 3;
- A pentapetide in water with Martini 2.3P polarizable force field[^5];
- Small molecules (Martini 3);
- Several _ad hoc_ systems to test specific potential terms.

All the test systems files are in the `tests` directory, which also includes example scripts to compare GROMACS and OpenMM energies and forces.

## Limitations
- Martini 2 cholesterol. The standard Martini 2 cholesterol topology uses a constraint network that, while can be solved by LINCS in GROMACS, cannot be solved by the Constant Constraint Matrix Approximation (CCMA) algorithm in OpenMM, thus leading to instabilities.
    - We have developed an alternative topology to be used in OpenMM. Before running an OpenMM Martini 2 simulation with cholesterol, the user needs to manually replace the cholesterol topology. The modified topology can be found in the folder `cholesterol`.
- Math in .itp files:
    - Gromacs allows for mathematical expressions to be used in .itp files. For example:
        - `    1     6      1   0.98112 RUBBER_FC*1.000000`
    - This is not allowed in OpenMM and must be edited:
        - `    1     6      1   0.98112 RUBBER_FC`

## Dependencies
- Tested with Python 3.8
- OpenMM v7.7.0

## Installation
To install the package, in your current Python environment run:

`python setup.py install`

## License
The source code included in this repository is available under the GNU Public License, version 3 (see `LICENSE`).

## Contributions
Please report any problems or bugs you may encounter as issues on [the github repo](https://github.com/maccallumlab/martini_openmm).

## References
[^1]: Marrink, S. J. et al., J. Phys. Chem. B 2007, 111, 7812–7824
[^2]: de Jong, D. H. et al., J. Chem. Theory Comput. 2013, 9, 687–697
[^3]: Souza, P. C. T. et al., Nat. Methods 2021, 18, 382–388
[^4]: Eastman, P. et al., PLoS Comput. Biol. 2017, 13, e1005659
[^5]: Khan, H. M. et al., J. Chem. Theory Comput. 2020, 16, 2550–2560
