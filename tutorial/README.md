# Tutorial: Using Martini force field in OpenMM simulations

*Shangnong Hu, University of Calgary, Calgary, Alberta, Canada*

## Introduction

This tutorial will walk you through how to set up an OpenMM[^1] simulation using the Martini force field[^2][^3][^4] and GROMACS[^5] topology files. More specifically, we use the [implementation](https://github.com/maccallumlab/martini_openmm) of Martini for OpenMM from the MacCallum Lab. If you are not familiar with the protocols for generating gro, top and itp files for Martini in GROMACS, please have a look in the tutorials from the [Martini official website](http://cgmartini.nl/) and [GROMACS manual](https://manual.gromacs.org/current/reference-manual/). Some familiarity with OpenMM is also assumed here. For explanation of how to work with OpenMM and run parameters we refer to the [OpenMM user guide](http://docs.openmm.org/latest/userguide/). There are some excellent OpenMM tutorials available at the same page as well. Although OpenMM supports Python and C++ as APIs, here we only use Python for demonstration.

We assume you are using a Linux environment. To install Python libraries for OpenMM and Martini, we recommend using a conda environment. You can install Miniconda or Anaconda following the instructions at [conda.io](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html). To create a new environment for Python development using conda, you can use the following command:

```bash
$ conda create --name myOpenMM-env python=3.8
```

Here we specifically request the Python version 3.8 because the Martini implementation is developed using that version. After creating a new environment, run the following commands to activate that environment and install OpenMM with all the dependencies and MDTraj[^6]:

```bash
$ conda activate myOpenMM-env
(myOpenMM-env)$ conda install -c conda-forge openmm
(myOpenMM-env)$ conda install -c conda-forge mdtraj
```

When an environment is successfully activated, the name of the environment will be displayed between parentheses at the left of each command line prompt. If you only see `(base)`, the environment is not correctly activated.

Then we will install the Martini implementation from the repo by cloning the `martini3` [branch](https://github.com/maccallumlab/martini_openmm/tree/martini3) to your local machine. In the conda environment you created earlier, run the following command to install the implementation:

```bash
(myOpenMM-env)$ python setup.py install
```

To follow this tutorial, you can download the sample script and files here. The tutorial folder includes all the files you need. We will show a brief workflow on how to perform Martini simulation in OpenMM. To explore further capabilities of OpenMM, you are invited to visit the official documentation.

## Setting up

In the example python script named run.py, we will simulate a ubiquitin in water. The structure was prepared with `martinize2`[^7] and GROMACS. At the beginning of the script, we import the following libraries:

```python
from openmm.unit import *
from openmm import *
from openmm.app import *
import martini_openmm as martini
from mdtraj.reporters import XTCReporter
from sys import stdout
```
In the beginning of the script, we specify whether we will use only CPU or mainly GPU. You also have the freedom to decide the computing precision. From line 10 to 31, we import the structure and topology in GROMACS format. We import the structure with a built-in GRO file reader from OpenMM:

```python
conf = GromacsGroFile("system.gro")
```

Then, we import the topology in Martini model:

```python
top = martini.MartiniTopFile("system.top",
		periodicBoxVectors=box_vectors,
		defines=defines,
		epsilon_r=epsilon_r
		)
```

Before minimizing the system, we define the non-bonded cut-off, integrator, target temperature, integration time steps. This information is stored in a Simulation object that has all the parameters.

## Minimization

The built-in function `minimizeEnergy()` is used along with the maximum iteration number and tolerance for energy convergence, which depend on the simulation system:

```python
simulation.minimizeEnergy(maxIterations=5000,tolerance=1.0)
```

## Equilibration

In our example, we perform a 1 ns NVT equilibration followed by a 1 ns NPT equilibration. From the output, you can notice the change of simulation box volume.

In NPT equilibration, we introduce a Monte-Carlo barostat for the pressure coupling:

```python
system.addForce(MonteCarloBarostat(1 * bar, 310 * kelvin))
```

At the end of the equilibration, we can save the simulation system into files:

```python
simulation.saveState('equi.state')
simulation.saveCheckpoint('equi.chk')
```

Knowing that PDB and other file formats are also available.

## Production

### Simulations in NPT ensemble

Here we use a Monte-Carlo barostat and a Langevin integrator to maintain pressure and temperature in our NPT simulation. We recommend using no larger than 20 fs time steps. The random number seed is set to 0 to have different forces at each run, or you can specify your own seed number across different simulations if you need to do so.

In the last section, we define the output files and the length of simulation. We use XTCReporter from mdtraj to generate the trajectory in XTC format:

```python
simulation.reporters.append(StateDataReporter("prod.log", 1000,
		step=True,
		potentialEnergy=True,
		totalEnergy=True,
		density=True,
		temperature=True,
		volume=True)
		)
xtc_reporter = XTCReporter('prod.xtc', 1000)
simulation.reporters.append(xtc_reporter)
```

Unlike GROMACS, the data wanted in the output file needs to be specified prior to the simulation. The list of properties available includes temperature, energies, density, volumes, etc…

#### Aqueous system

After comparison with GROMACS simulation, the friction coefficient in Langevin integrator in OpenMM produces the closest result when set to 10 per picosecond.

#### Lipid bilayer system

When simulating a lipid bilayer system, the barostat needs to be decoupled between the bilayer plane (often X and Y axis) and the vertical axis (often Z). In OpenMM, we need to use a dedicated Monte-Carlo membrane barostat:

```python
barostat = mm.openmm.MonteCarloMembraneBarostat(1 * bar, 0 * bar * nanometer, 310 * kelvin,
		mm.openmm.MonteCarloMembraneBarostat.XYIsotropic,
		mm.openmm.MonteCarloMembraneBarostat.ZFree, 10
		)
```

### Simulations in NVT ensemble

To run an NVT simulation, we can use the same script as for NPT but omit the barostat. The rest remains the same.

### Simulations in NVE ensemble

Although NVE simulation has less frequent usage than the previous two ensembles, we will mention briefly its setup. In principle, we omit the barostat and the thermostat in our simulation script. In OpenMM, we will need to use the Velocity Verlet integrator by importing from openmmtools. You might need to install it before the first import. Then add the following line in the import section:

```python
from openmmtools.integrators import VelocityVerletIntegrator
```

By replacing the integrator with the new one, only the time steps need to be specified:

```python
integrator = VelocityVerletIntegrator(20 * femtosecond)
```

## Run

The sample script provided for this tutorial can be ran using the following command in the conda environment you created:

```bash
(myOpenMM-env)$ python run.py
```

## Analysis

The log file stores all the properties we specified along the simulation in a tabulated format. You can easily plot it using a data visualization software or parsing the file:

```
#"Step","Potential Energy (kJ/mole)","Total Energy (kJ/mole)","Temperature (K)"
200,-5164.372920680828,-3594.344823877926,328.6875138132883
400,-5138.979121536883,-3594.1487321152536,323.4123395601584
600,-5161.453751389345,-3594.844274865674,327.97182101610815
800,-5160.404523900361,-3595.146624271858,327.68886655796507
1000,-5159.04122180917,-3594.25071109566,327.59101805400155
1200,-5247.469255940026,-3595.130230175832,345.9194185511981
…
```

As the trajectory output is in XTC format, many of analysis proposed in the Martini tutorials using GROMACS can be used here as well. These tutorials are available in the official Martini website.

## References

[^1]:Eastman P, Swails J, Chodera JD, McGibbon RT, Zhao Y, Beauchamp KA, et al. OpenMM 7: Rapid development of high performance algorithms for molecular dynamics. PLOS Comput Biol. 2017 Jul 26;13(7):e1005659.
[^2]:Marrink SJ, Risselada HJ, Yefimov S, Tieleman DP, de Vries AH. The MARTINI Force Field: Coarse Grained Model for Biomolecular Simulations. J Phys Chem B. 2007 Jul 1;111(27):7812–24.
[^3]:Monticelli L, Kandasamy SK, Periole X, Larson RG, Tieleman DP, Marrink SJ. The MARTINI Coarse-Grained Force Field: Extension to Proteins. J Chem Theory Comput. 2008 May 1;4(5):819–34.
[^4]:Souza PCT, Alessandri R, Barnoud J, Thallmair S, Faustino I, Grünewald F, et al. Martini 3: a general purpose force field for coarse-grained molecular dynamics. Nat Methods. 2021 Apr;18(4):382–8.
[^5]:Abraham MJ, Murtola T, Schulz R, Páll S, Smith JC, Hess B, et al. GROMACS: High performance molecular simulations through multi-level parallelism from laptops to supercomputers. SoftwareX. 2015 Sep 1;1–2:19–25.
[^6]:McGibbon RT, Beauchamp KA, Harrigan MP, Klein C, Swails JM, Hernández CX, et al. MDTraj: A Modern Open Library for the Analysis of Molecular Dynamics Trajectories. Biophys J. 2015 Oct 20;109(8):1528–32.
[^7]:Kroon PC, Grünewald F, Barnoud J, van Tilburg M, Souza PCT, Wassenaar TA, et al. Martinize2 and Vermouth: Unified Framework for Topology Generation [Internet]. arXiv; 2022 [cited 2023 Feb 8]. Available from: http://arxiv.org/abs/2212.01191