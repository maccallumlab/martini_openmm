"""
martini.py: Used for loading Gromacs martini top files.

This is part of the OpenMM molecular simulation toolkit originating from
Simbios, the NIH National Center for Physics-Based Simulation of
Biological Structures at Stanford, funded under the NIH Roadmap for
Medical Research, grant U54 GM072970. See https://simtk.org.

Portions copyright (c) 2012-2018 Stanford University and the Authors.
Authors: Peter Eastman
Contributors: Jason Swails
Contributors: Justin MacCallum

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the "Software"),
to deal in the Software without restriction, including without limitation
the rights to use, copy, modify, merge, publish, distribute, sublicense,
and/or sell copies of the Software, and to permit persons to whom the
Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
USE OR OTHER DEALINGS IN THE SOFTWARE.
"""
__author__ = "Justin MacCallum"
__version__ = "0.1"

import distutils.spawn
import math
import os
import re
from collections import OrderedDict, defaultdict

import simtk.openmm as mm
import simtk.unit as unit
from simtk.openmm.app import PDBFile, Topology
from simtk.openmm.app import amberprmtopfile as prmtop
from simtk.openmm.app import element as elem
from simtk.openmm.app import forcefield as ff
from .vsites import (
    LinearSite,
    OutOfPlane,
    VSiteManager,
    COMLinearSite,
    NormalizedInPlaneSite,
    NormalizedInPlaneTwoParticleSite,
)

HBonds = ff.HBonds
AllBonds = ff.AllBonds
HAngles = ff.HAngles

OBC2 = prmtop.OBC2

novarcharre = re.compile(r"\W")


def _find_all_instances_in_string(string, substr):
    """Find indices of all instances of substr in string"""
    indices = []
    idx = string.find(substr, 0)
    while idx > -1:
        indices.append(idx)
        idx = string.find(substr, idx + 1)
    return indices


def _replace_defines(line, defines):
    """Replaces defined tokens in a given line"""
    if not defines:
        return line
    for define in reversed(defines):
        value = defines[define]
        indices = _find_all_instances_in_string(line, define)
        if not indices:
            continue
        # Check to see if it's inside of quotes
        inside = ""
        idx = 0
        n_to_skip = 0
        new_line = []
        for i, char in enumerate(line):
            if n_to_skip:
                n_to_skip -= 1
                continue
            if char in ("'\""):
                if not inside:
                    inside = char
                else:
                    if inside == char:
                        inside = ""
            if idx < len(indices) and i == indices[idx]:
                if inside:
                    new_line.append(char)
                    idx += 1
                    continue
                if i == 0 or novarcharre.match(line[i - 1]):
                    endidx = indices[idx] + len(define)
                    if endidx >= len(line) or novarcharre.match(line[endidx]):
                        new_line.extend(list(value))
                        n_to_skip = len(define) - 1
                        idx += 1
                        continue
                idx += 1
            new_line.append(char)
        line = "".join(new_line)

    return line


class MartiniTopFile(object):
    """Parse a Gromacs Martini top file and constructs a Topology and (optionally) an OpenMM System from it."""

    class _MoleculeType(object):
        """Inner class to store information about a molecule type."""

        def __init__(self):
            self.molecule_name = ""
            self.atoms = []
            self.bonds = []
            self.g96_angles = []
            self.harmonic_angles = []
            self.restricted_angles = []
            self.dihedrals = []
            self.exclusions = []
            self.pairs = []
            self.constraints = []
            self.cmaps = []
            self.vsites = VSiteManager()

    def __init__(
        self,
        file,
        periodicBoxVectors=None,
        unitCellDimensions=None,
        includeDir=None,
        defines=None,
        epsilon_r=15.0,
    ):
        """Load a top file.

        Parameters
        ----------
        file : str
            the name of the file to load
        periodicBoxVectors : tuple of Vec3=None
            the vectors defining the periodic box
        unitCellDimensions : Vec3=None
            the dimensions of the crystallographic unit cell.  For
            non-rectangular unit cells, specify periodicBoxVectors instead.
        includeDir : string=None
            A directory in which to look for other files included from the
            top file. If not specified, we will attempt to locate a gromacs
            installation on your system. When gromacs is installed in
            /usr/local, this will resolve to /usr/local/gromacs/share/gromacs/top
        defines : dict={}
            preprocessor definitions that should be predefined when parsing the file
        """
        self.epsilon_r = epsilon_r
        if includeDir is None:
            includeDir = _get_default_gromacs_include_dir()
        self._includeDirs = (os.path.dirname(file), includeDir)
        self._defines = OrderedDict()
        self._genpairs = True
        if defines is not None:
            for define, value in defines.items():
                self._defines[define] = value

        self._currentCategory = None
        self._ifStack = []
        self._elseStack = []
        self._moleculeTypes = {}
        self._molecules = []
        self._currentMoleculeType = None
        self._atom_types = {}
        self._bondTypes = {}
        self._angleTypes = {}
        self._dihedralTypes = {}
        self._implicitTypes = {}
        self._pairTypes = {}
        self._cmapTypes = {}
        self._nonbond_types = {}
        self._all_vsites = []
        self.nb_force = None
        self.es_self_excl_force = None
        self.es_except_force = None
        self.lj_except_force = None
        self.harmonic_angle_force = None
        self.g96_angle_force = None
        self.restricted_angle_force = None
        self.harmonic_bond_force = None
        self.periodic_torsion_force = None
        self.harmonic_torsion_force = None
        self.combined_bending_torsion_force = None
        self.rb_torsion_force = None
        self._use_sigma_eps = False
        self._use_g96_angles = False
        self._use_harmonic_angles = False
        self._use_restricted_angles = False

        self._process_file(file)

        top = Topology()
        self.topology = top
        if periodicBoxVectors is not None:
            if unitCellDimensions is not None:
                raise ValueError(
                    "specify either periodicBoxVectors or unitCellDimensions, but not both"
                )
            top.setPeriodicBoxVectors(periodicBoxVectors)
        else:
            top.setUnitCellDimensions(unitCellDimensions)
        PDBFile._loadNameReplacementTables()
        for molecule_name, molecule_count in self._molecules:
            if molecule_name not in self._moleculeTypes:
                raise ValueError("Unknown molecule type: " + molecule_name)
            molecule_type = self._moleculeTypes[molecule_name]

            # Create the specified number of molecules of this type.
            for i in range(molecule_count):
                atoms = []
                lastResidue = None
                c = top.addChain()
                for index, fields in enumerate(molecule_type.atoms):
                    resNumber = fields[2]
                    if resNumber != lastResidue:
                        lastResidue = resNumber
                        resName = fields[3]
                        if resName in PDBFile._residueNameReplacements:
                            resName = PDBFile._residueNameReplacements[resName]
                        r = top.addResidue(resName, c)
                        if resName in PDBFile._atomNameReplacements:
                            atomReplacements = PDBFile._atomNameReplacements[resName]
                        else:
                            atomReplacements = {}
                    atomName = fields[4]
                    if atomName in atomReplacements:
                        atomName = atomReplacements[atomName]

                    # Try to guess the element.
                    upper = atomName.upper()
                    if upper.startswith("CL"):
                        element = elem.chlorine
                    elif upper.startswith("NA"):
                        element = elem.sodium
                    elif upper.startswith("MG"):
                        element = elem.magnesium
                    else:
                        try:
                            element = elem.get_by_symbol(atomName[0])
                        except KeyError:
                            element = None
                    atoms.append(top.addAtom(atomName, element, r))

                # Add bonds to the topology

                for fields in molecule_type.bonds:
                    top.addBond(atoms[int(fields[0]) - 1], atoms[int(fields[1]) - 1])

    def _process_file(self, file):
        append = ""
        with open(file) as lines:
            for line in lines:
                if line.strip().endswith("\\"):
                    append = "%s %s" % (append, line[: line.rfind("\\")])
                else:
                    self._process_line(append + " " + line, file)
                    append = ""

    def _process_line(self, line, file):
        """Process one line from a file."""
        if ";" in line:
            line = line[: line.index(";")]
        stripped = line.strip()
        ignore = not all(self._ifStack)
        if stripped.startswith("*") or len(stripped) == 0:
            # A comment or empty line.
            return

        elif stripped.startswith("[") and not ignore:
            # The start of a category.
            if not stripped.endswith("]"):
                raise ValueError("Illegal line in .top file: " + line)
            self._currentCategory = stripped[1:-1].strip()

        elif stripped.startswith("#"):
            # A preprocessor command.
            fields = stripped.split()
            command = fields[0]
            if len(self._ifStack) != len(self._elseStack):
                raise RuntimeError("#if/#else stack out of sync")

            if command == "#include" and not ignore:
                # Locate the file to include
                name = stripped[len(command) :].strip(' \t"<>')
                searchDirs = self._includeDirs + (os.path.dirname(file),)
                for dir in searchDirs:
                    file = os.path.join(dir, name)
                    if os.path.isfile(file):
                        # We found the file, so process it.
                        self._process_file(file)
                        break
                else:
                    raise ValueError("Could not locate #include file: " + name)
            elif command == "#define" and not ignore:
                # Add a value to our list of defines.
                if len(fields) < 2:
                    raise ValueError("Illegal line in .top file: " + line)
                name = fields[1]
                valueStart = stripped.find(name, len(command)) + len(name) + 1
                value = line[valueStart:].strip()
                value = value or "1"  # Default define is 1
                self._defines[name] = value
            elif command == "#ifdef":
                # See whether this block should be ignored.
                if len(fields) < 2:
                    raise ValueError("Illegal line in .top file: " + line)
                name = fields[1]
                self._ifStack.append(name in self._defines)
                self._elseStack.append(False)
            elif command == "#undef":
                # Un-define a variable
                if len(fields) < 2:
                    raise ValueError("Illegal line in .top file: " + line)
                if fields[1] in self._defines:
                    self._defines.pop(fields[1])
            elif command == "#ifndef":
                # See whether this block should be ignored.
                if len(fields) < 2:
                    raise ValueError("Illegal line in .top file: " + line)
                name = fields[1]
                self._ifStack.append(name not in self._defines)
                self._elseStack.append(False)
            elif command == "#endif":
                # Pop an entry off the if stack.
                if len(self._ifStack) == 0:
                    raise ValueError("Unexpected line in .top file: " + line)
                del self._ifStack[-1]
                del self._elseStack[-1]
            elif command == "#else":
                # Reverse the last entry on the if stack
                if len(self._ifStack) == 0:
                    raise ValueError("Unexpected line in .top file: " + line)
                if self._elseStack[-1]:
                    raise ValueError(
                        "Unexpected line in .top file: "
                        "#else has already been used " + line
                    )
                self._ifStack[-1] = not self._ifStack[-1]
                self._elseStack[-1] = True

        elif not ignore:
            # Gromacs occasionally uses #define's to introduce specific
            # parameters for individual terms (for instance, this is how
            # ff99SB-ILDN is implemented). So make sure we do the appropriate
            # pre-processor replacements necessary
            line = _replace_defines(line, self._defines)
            # A line of data for the current category
            if self._currentCategory is None:
                raise ValueError("Unexpected line in .top file: " + line)
            if self._currentCategory == "defaults":
                self._process_defaults(line)
            elif self._currentCategory == "moleculetype":
                self._process_molecule_type(line)
            elif self._currentCategory == "molecules":
                self._process_molecule(line)
            elif self._currentCategory == "atoms":
                self._process_atom(line)
            elif self._currentCategory == "bonds":
                self._process_bond(line)
            elif self._currentCategory == "angles":
                self._process_angle(line)
            elif self._currentCategory == "dihedrals":
                self._process_dihedral(line)
            elif self._currentCategory == "exclusions":
                self._process_exclusion(line)
            elif self._currentCategory == "pairs":
                self._process_pair(line)
            elif self._currentCategory == "constraints":
                self._process_constraint(line)
            elif self._currentCategory == "cmap":
                self._process_cmap(line)
            elif self._currentCategory == "atomtypes":
                self._process_atom_type(line)
            elif self._currentCategory == "bondtypes":
                self._process_bond_type(line)
            elif self._currentCategory == "angletypes":
                self._process_angle_type(line)
            elif self._currentCategory == "dihedraltypes":
                self._process_dihedral_type(line)
            elif self._currentCategory == "implicit_genborn_params":
                self._process_implicit_type(line)
            elif self._currentCategory == "pairtypes":
                self._process_pair_type(line)
            elif self._currentCategory == "cmaptypes":
                self._process_cmap_type(line)
            elif self._currentCategory == "nonbond_params":
                self._process_nonbond_type(line)
            elif self._currentCategory == "virtual_sites2":
                self._process_vsites2(line)
            elif self._currentCategory == "virtual_sites3":
                self._process_vsites3(line)
            elif self._currentCategory == "virtual_sitesn":
                self._process_vsitesn(line)
            elif self._currentCategory.startswith("virtual_sites"):
                if self._currentMoleculeType is None:
                    raise ValueError(
                        "Found %s before [ moleculetype ]" % self._currentCategory
                    )
                raise ValueError(
                    f"[ {self._currentCategory} ] is currently not supported.\n"
                    "Please file an issue at github.com/maccallumlab/martini_openmm"
                )
            elif self._currentCategory == "system":  # ignore the system description
                pass
            elif self._currentCategory == "position_restraints":
                raise ValueError(
                    "[ position_restraints ] is not currently supported.\n"
                    "The same effect can be achieved using a CustomExternalForce in OpenMM.\n"
                    "Please see github.com/maccallumlab/martini_openmm for details."
                )
            else:
                raise ValueError(
                    f"[ {self._currentCategory} ] is not currently supported.\n"
                    "Please file an issue at github.com/maccallumlab/martini_openmm"
                )

    def _process_defaults(self, line):
        """Process the [ defaults ] line."""

        fields = line.split()
        if len(fields) > 2:
            raise ValueError("Too many fields in [ defaults ] line: " + line)
        if fields[0] != "1":
            raise ValueError("Unsupported nonbonded type: " + fields[0])
        if fields[1] == "1":
            self._use_sigma_eps = False
        elif fields[1] == "2":
            self._use_sigma_eps = True
        else:
            raise ValueError("Unsupported combination rule: " + fields[1])
        self._defaults = fields

    def _process_molecule_type(self, line):
        """Process a line in the [ moleculetypes ] category."""
        fields = line.split()
        if len(fields) < 1:
            raise ValueError("Too few fields in [ moleculetypes ] line: " + line)
        type = MartiniTopFile._MoleculeType()
        type.molecule_name = fields[0]
        self._moleculeTypes[fields[0]] = type
        self._currentMoleculeType = type

    def _process_molecule(self, line):
        """Process a line in the [ molecules ] category."""
        fields = line.split()
        if len(fields) < 2:
            raise ValueError("Too few fields in [ molecules ] line: " + line)
        self._molecules.append((fields[0], int(fields[1])))

    def _process_atom(self, line):
        """Process a line in the [ atoms ] category."""
        if self._currentMoleculeType is None:
            raise ValueError("Found [ atoms ] section before [ moleculetype ]")
        fields = line.split()
        if len(fields) < 5:
            raise ValueError("Too few fields in [ atoms ] line: " + line)
        self._currentMoleculeType.atoms.append(fields)

    def _process_bond(self, line):
        """Process a line in the [ bonds ] category."""
        if self._currentMoleculeType is None:
            raise ValueError("Found [ bonds ] section before [ moleculetype ]")
        fields = line.split()
        if len(fields) < 3:
            raise ValueError("Too few fields in [ bonds ] line: " + line)
        if fields[2] != "1" and fields[2] != "6":
            raise ValueError("Unsupported function type in [ bonds ] line: " + line)
        self._currentMoleculeType.bonds.append(fields)

    def _process_angle(self, line):
        """Process a line in the [ angles ] category."""
        if self._currentMoleculeType is None:
            raise ValueError("Found [ angles ] section before [ moleculetype ]")
        fields = line.split()

        if len(fields) < 4:
            raise ValueError("Too few fields in [ angles ] line: " + line)

        if fields[3] == "1":
            self._use_harmonic_angles = True
            self._currentMoleculeType.harmonic_angles.append(fields)
        elif fields[3] == "2":
            self._use_g96_angles = True
            self._currentMoleculeType.g96_angles.append(fields)
        elif fields[3] == "10":
            self._use_restricted_angles = True
            self._currentMoleculeType.restricted_angles.append(fields)
        else:
            raise ValueError(
                "Unsupported (nonG96) function type in [ angles ] line: " + line
            )

    def _process_dihedral(self, line):
        """Process a line in the [ dihedrals ] category."""
        if self._currentMoleculeType is None:
            raise ValueError("Found [ dihedrals ] section before [ moleculetype ]")
        fields = line.split()
        if len(fields) < 5:
            raise ValueError("Too few fields in [ dihedrals ] line: " + line)
        if fields[4] not in ("1", "2", "3", "4", "5", "9", "11"):
            raise ValueError("Unsupported function type in [ dihedrals ] line: " + line)
        self._currentMoleculeType.dihedrals.append(fields)

    def _process_exclusion(self, line):
        """Process a line in the [ exclusions ] category."""
        if self._currentMoleculeType is None:
            raise ValueError("Found [ exclusions ] section before [ moleculetype ]")
        fields = line.split()
        if len(fields) < 2:
            raise ValueError("Too few fields in [ exclusions ] line: " + line)
        self._currentMoleculeType.exclusions.append(fields)

    def _process_pair(self, line):
        """Process a line in the [ pairs ] category."""
        if self._currentMoleculeType is None:
            raise ValueError("Found [ pairs ] section before [ moleculetype ]")
        fields = line.split()
        if len(fields) < 3:
            raise ValueError("Too few fields in [ pairs ] line: " + line)
        if fields[2] != "1":
            raise ValueError("Unsupported function type in [ pairs ] line: " + line)
        self._currentMoleculeType.pairs.append(fields)

    def _process_constraint(self, line):
        """Process a line in the [ constraints ] category."""
        if self._currentMoleculeType is None:
            raise ValueError("Found [ constraints ] section before [ moleculetype ]")
        fields = line.split()
        if len(fields) < 4:
            raise ValueError("Too few fields in [ constraints ] line: " + line)
        self._currentMoleculeType.constraints.append(fields)

    def _process_cmap(self, line):
        """Process a line in the [ cmaps ] category."""
        if self._currentMoleculeType is None:
            raise ValueError("Found [ cmap ] section before [ moleculetype ]")
        fields = line.split()
        if len(fields) < 6:
            raise ValueError("Too few fields in [ cmap ] line: " + line)
        self._currentMoleculeType.cmaps.append(fields)

    def _process_atom_type(self, line):
        """Process a line in the [ atomtypes ] category."""
        fields = line.split()
        if len(fields) < 6:
            raise ValueError("Too few fields in [ atomtypes ] line: " + line)
        if len(fields[3]) == 1:
            # Bonded type and atomic number are both missing.
            fields.insert(1, None)
            fields.insert(1, None)
        elif len(fields[4]) == 1 and fields[4].isalpha():
            if fields[1][0].isalpha():
                # Atomic number is missing.
                fields.insert(2, None)
            else:
                # Bonded type is missing.
                fields.insert(1, None)
        self._atom_types[fields[0]] = fields

    def _process_bond_type(self, line):
        """Process a line in the [ bondtypes ] category."""
        fields = line.split()
        if len(fields) < 5:
            raise ValueError("Too few fields in [ bondtypes ] line: " + line)
        if fields[2] != "1":
            raise ValueError("Unsupported function type in [ bondtypes ] line: " + line)
        self._bondTypes[tuple(fields[:2])] = fields

    def _process_angle_type(self, line):
        """Process a line in the [ angletypes ] category."""
        fields = line.split()
        if len(fields) < 6:
            raise ValueError("Too few fields in [ angletypes ] line: " + line)
        if fields[3] not in ("1", "5"):
            raise ValueError(
                "Unsupported function type in [ angletypes ] line: " + line
            )
        self._angleTypes[tuple(fields[:3])] = fields

    def _process_dihedral_type(self, line):
        """Process a line in the [ dihedraltypes ] category."""
        fields = line.split()
        if len(fields) < 7:
            raise ValueError("Too few fields in [ dihedraltypes ] line: " + line)
        if fields[4] not in ("1", "2", "3", "4", "5", "9"):
            raise ValueError(
                "Unsupported function type in [ dihedraltypes ] line: " + line
            )
        key = tuple(fields[:5])
        if fields[4] == "9" and key in self._dihedralTypes:
            # There are multiple dihedrals defined for these atom types.
            self._dihedralTypes[key].append(fields)
        else:
            self._dihedralTypes[key] = [fields]

    def _process_implicit_type(self, line):
        """Process a line in the [ implicit_genborn_params ] category."""
        fields = line.split()
        if len(fields) < 6:
            raise ValueError(
                "Too few fields in [ implicit_genborn_params ] line: " + line
            )
        self._implicitTypes[fields[0]] = fields

    def _process_pair_type(self, line):
        """Process a line in the [ pairtypes ] category."""
        fields = line.split()
        if len(fields) < 5:
            raise ValueError("Too few fields in [ pairtypes] line: " + line)
        if fields[2] != "1":
            raise ValueError("Unsupported function type in [ pairtypes ] line: " + line)
        self._pairTypes[tuple(fields[:2])] = fields

    def _process_cmap_type(self, line):
        """Process a line in the [ cmaptypes ] category."""
        fields = line.split()
        if len(fields) < 8 or len(fields) < 8 + int(fields[6]) * int(fields[7]):
            raise ValueError("Too few fields in [ cmaptypes ] line: " + line)
        if fields[5] != "1":
            raise ValueError("Unsupported function type in [ cmaptypes ] line: " + line)
        self._cmapTypes[tuple(fields[:5])] = fields

    def _process_nonbond_type(self, line):
        """Process a line in the [ nonbond_params ] category."""
        fields = line.split()
        if len(fields) < 5:
            raise ValueError("Too few fields in [ nonbond_params ] line: " + line)
        if fields[2] != "1":
            raise ValueError(
                "Unsupported function type in [ nonbond_params ] line: " + line
            )
        self._nonbond_types[tuple(sorted(fields[:2]))] = fields

    def _process_vsites2(self, line):
        """Process a line in the [ virtual_sites2 ] category."""
        fields = line.split()
        if len(fields) < 5:
            raise ValueError(f"Too few fields in [ virtual_sites2 ] line: {line}")

        index = int(fields[0])
        atom1 = int(fields[1])
        atom2 = int(fields[2])
        func = int(fields[3])
        if func == 1:
            w = float(fields[4])

            site_dict = {atom1: 1 - w, atom2: w}
            site = LinearSite(site_dict)

            self._currentMoleculeType.vsites.add(index, site)
        elif func == 2:
            a = float(fields[4])
            site = NormalizedInPlaneTwoParticleSite(atom1, atom2, a)
            self._currentMoleculeType.vsites.add(index, site)
        else:
            raise ValueError(
                f"Unknown function type {func} in virtual_sites2 line: {line}"
            )

    def _process_vsites3(self, line):
        """Process a line in the [ virtual_sites3 ] category."""
        fields = line.split()
        if len(fields) < 7:
            raise ValueError(f"Too few fields in [ virtual_sites3 ] line: {line}")

        index = int(fields[0])
        atom1 = int(fields[1])
        atom2 = int(fields[2])
        atom3 = int(fields[3])
        func = int(fields[4])

        # type 3 in gromacs
        if func == 1:
            if len(fields) < 7:
                raise ValueError(f"Not enough parameters for vsite: {line}")
            a = float(fields[5])
            b = float(fields[6])
            w1 = 1 - a - b
            w2 = a
            w3 = b
            site_dict = {
                atom1: w1,
                atom2: w2,
                atom3: w3,
            }
            site = LinearSite(site_dict)
            self._currentMoleculeType.vsites.add(index, site)
        # type 3fd in gromacs
        elif func == 2:
            if len(fields) < 7:
                raise ValueError(f"Not enough parameters for vsite: {line}")
            a = float(fields[5])
            d = float(fields[6])
            site = NormalizedInPlaneSite(atom1, atom2, atom3, a, d)
            self._currentMoleculeType.vsites.add(index, site)
        # type 3out in gromacs
        elif func == 4:
            if len(fields) < 8:
                raise ValueError(f"Not enough parameters for type 4 site: {line}")
            a = float(fields[5])
            b = float(fields[6])
            c = float(fields[7])
            site = OutOfPlane(atom1, atom2, atom3, a, b, c)
            self._currentMoleculeType.vsites.add(index, site)
        else:
            raise RuntimeError(f"Function type {func} unsupported in vsites3: {line}")

    def _process_vsitesn(self, line):
        """Process a line in the [ virtual_sitesn ] category."""
        fields = line.split()
        index = int(fields[0])
        func = int(fields[1])
        from_atoms = [int(field) for field in fields[2:]]

        if len(from_atoms) == 1:
            site_dict = {from_atoms[0]: 1.0}
            site = LinearSite(site_dict)
            self._currentMoleculeType.vsites.add(index, site)

        elif func == 1:  # center of geometry
            w = 1 / len(from_atoms)
            site_dict = {atom: w for atom in from_atoms}
            site = LinearSite(site_dict)
            self._currentMoleculeType.vsites.add(index, site)

        elif func == 2:  # center of mass
            site = COMLinearSite(from_atoms)
            self._currentMoleculeType.vsites.add(index, site)

        else:
            raise ValueError(f"Site type {func} unspported: {line}")

    def _build_dihedral_lookup_table(self):
        dihedral_type_table = {}
        for key in self._dihedralTypes:
            if key[1] != "X" and key[2] != "X":
                if (key[1], key[2]) not in dihedral_type_table:
                    dihedral_type_table[(key[1], key[2])] = []
                dihedral_type_table[(key[1], key[2])].append(key)
                if (key[2], key[1]) not in dihedral_type_table:
                    dihedral_type_table[(key[2], key[1])] = []
                dihedral_type_table[(key[2], key[1])].append(key)
        wildcard_dihedral_types = []
        for key in self._dihedralTypes:
            if key[1] == "X" or key[2] == "X":
                wildcard_dihedral_types.append(key)
                for types in dihedral_type_table.values():
                    types.append(key)
        return dihedral_type_table, wildcard_dihedral_types

    def _add_atoms_to_system(self, sys, molecule_type):
        for fields in molecule_type.atoms:
            if len(fields) >= 8:
                mass = float(fields[7])
            else:
                mass = float(self._atom_types[fields[1]][3])
            sys.addParticle(mass)

    def _add_bonds_to_system(self, sys, molecule_type, bonded_types, base_atom_index):
        for fields in molecule_type.bonds:
            atoms = [int(x) - 1 for x in fields[:2]]
            types = tuple(bonded_types[i] for i in atoms)
            if len(fields) >= 5:
                params = fields[3:5]
            elif types in self._bondTypes:
                params = self._bondTypes[types][3:5]
            elif types[::-1] in self._bondTypes:
                params = self._bondTypes[types[::-1]][3:5]
            else:
                raise ValueError(
                    "No parameters specified for bond: " + fields[0] + ", " + fields[1]
                )

            length = float(params[0])

            if self.harmonic_bond_force is None:
                self.harmonic_bond_force = mm.HarmonicBondForce()
                sys.addForce(self.harmonic_bond_force)
            self.harmonic_bond_force.addBond(
                base_atom_index + atoms[0],
                base_atom_index + atoms[1],
                length,
                float(params[1]),
            )

    def _add_angles_to_system(self, sys, molecule_type, bonded_types, base_atom_index):
        self._add_angle_forces(sys)
        self._add_harmonic_angles_to_system(
            molecule_type, bonded_types, base_atom_index
        )
        self._add_g96_angles_to_system(molecule_type, bonded_types, base_atom_index)
        self._add_restricted_angles_to_system(
            molecule_type, bonded_types, base_atom_index
        )

    def _add_angle_forces(self, sys):
        if self.harmonic_angle_force is None:
            self.harmonic_angle_force = mm.HarmonicAngleForce()
            sys.addForce(self.harmonic_angle_force)

        if self._use_g96_angles and self.g96_angle_force is None:
            self.g96_angle_force = mm.CustomAngleForce(
                "0.5 * k * (cos(theta) - cos(theta0))^2"
            )
            self.g96_angle_force.addPerAngleParameter("theta0")
            self.g96_angle_force.addPerAngleParameter("k")
            sys.addForce(self.g96_angle_force)

        if self._use_restricted_angles and self.restricted_angle_force is None:
            self.restricted_angle_force = mm.CustomAngleForce(
                "0.5 * k * (cos(theta) - cos(theta0))^2 / sin(theta)^2"
            )
            self.restricted_angle_force.addPerAngleParameter("theta0")
            self.restricted_angle_force.addPerAngleParameter("k")
            sys.addForce(self.restricted_angle_force)

    def _process_angle_params(self, fields, bonded_types):
        degToRad = math.pi / 180
        atoms = [int(x) - 1 for x in fields[:3]]
        types = tuple(bonded_types[i] for i in atoms)
        if len(fields) >= 6:
            params = fields[4:]
        elif types in self._angleTypes:
            params = self._angleTypes[types][4:]
        elif types[::-1] in self._angleTypes:
            params = self._angleTypes[types[::-1]][4:]
        else:
            raise ValueError(
                "No parameters specified for angle: "
                + fields[0]
                + ", "
                + fields[1]
                + ", "
                + fields[2]
            )

        theta = float(params[0]) * degToRad
        return atoms, theta, float(params[1])

    def _add_harmonic_angles_to_system(
        self, molecule_type, bonded_types, base_atom_index
    ):
        for fields in molecule_type.harmonic_angles:
            assert fields[3] == "1"
            atoms, theta, k = self._process_angle_params(fields, bonded_types)
            self.harmonic_angle_force.addAngle(
                base_atom_index + atoms[0],
                base_atom_index + atoms[1],
                base_atom_index + atoms[2],
                theta,
                k,
            )

    def _add_g96_angles_to_system(self, molecule_type, bonded_types, base_atom_index):
        for fields in molecule_type.g96_angles:
            assert fields[3] == "2"
            atoms, theta, k = self._process_angle_params(fields, bonded_types)
            self.g96_angle_force.addAngle(
                base_atom_index + atoms[0],
                base_atom_index + atoms[1],
                base_atom_index + atoms[2],
                [theta, k],
            )

    def _add_restricted_angles_to_system(
        self, molecule_type, bonded_types, base_atom_index
    ):
        for fields in molecule_type.restricted_angles:
            assert fields[3] == "10"
            atoms, theta, k = self._process_angle_params(fields, bonded_types)
            self.restricted_angle_force.addAngle(
                base_atom_index + atoms[0],
                base_atom_index + atoms[1],
                base_atom_index + atoms[2],
                [theta, k],
            )

    def _add_torsions_to_system(
        self,
        sys,
        molecule_type,
        bonded_types,
        dihedral_type_table,
        wildcard_dihedral_types,
        base_atom_index,
    ):
        degToRad = math.pi / 180

        for fields in molecule_type.dihedrals:
            atoms = [int(x) - 1 for x in fields[:4]]
            types = tuple(bonded_types[i] for i in atoms)
            dihedralType = fields[4]
            reversedTypes = types[::-1] + (dihedralType,)
            types = types + (dihedralType,)
            if (
                (dihedralType in ("1", "4", "5", "9", "11") and len(fields) > 7)
                or (dihedralType == "3" and len(fields) > 10)
                or (dihedralType == "2" and len(fields) > 6)
            ):
                paramsList = [fields]
            else:
                # Look for a matching dihedral type.
                paramsList = None
                if (types[1], types[2]) in dihedral_type_table:
                    dihedralTypes = dihedral_type_table[(types[1], types[2])]
                else:
                    dihedralTypes = wildcard_dihedral_types
                for key in dihedralTypes:
                    if all(a == b or a == "X" for a, b in zip(key, types)) or all(
                        a == b or a == "X" for a, b in zip(key, reversedTypes)
                    ):
                        paramsList = self._dihedralTypes[key]
                        if "X" not in key:
                            break
                if paramsList is None:
                    raise ValueError(
                        "No parameters specified for dihedral: "
                        + fields[0]
                        + ", "
                        + fields[1]
                        + ", "
                        + fields[2]
                        + ", "
                        + fields[3]
                    )
            for params in paramsList:
                if dihedralType in ("1", "4", "9"):
                    # Periodic torsion
                    k = float(params[6])
                    if k != 0:
                        if self.periodic_torsion_force is None:
                            self.periodic_torsion_force = mm.PeriodicTorsionForce()
                            sys.addForce(self.periodic_torsion_force)
                        self.periodic_torsion_force.addTorsion(
                            base_atom_index + atoms[0],
                            base_atom_index + atoms[1],
                            base_atom_index + atoms[2],
                            base_atom_index + atoms[3],
                            int(float(params[7])),
                            float(params[5]) * degToRad,
                            k,
                        )
                elif dihedralType == "2":
                    # Harmonic torsion
                    k = float(params[6])
                    phi0 = float(params[5])
                    if k != 0:
                        if self.harmonic_torsion_force is None:
                            self.harmonic_torsion_force = mm.CustomTorsionForce(
                                "0.5*k*(thetap-theta0)^2; thetap = step(-(theta-theta0+pi))*2*pi+theta+step(theta-theta0-pi)*(-2*pi); pi = %.15g"
                                % math.pi
                            )
                            self.harmonic_torsion_force.addPerTorsionParameter("theta0")
                            self.harmonic_torsion_force.addPerTorsionParameter("k")
                            sys.addForce(self.harmonic_torsion_force)
                        # map phi0 into correct space
                        phi0 = phi0 - 360 if phi0 > 180 else phi0
                        self.harmonic_torsion_force.addTorsion(
                            base_atom_index + atoms[0],
                            base_atom_index + atoms[1],
                            base_atom_index + atoms[2],
                            base_atom_index + atoms[3],
                            (phi0 * degToRad, k),
                        )
                elif dihedralType == "11":
                    # combined bending / torsion
                    k = float(params[5])
                    a0 = float(params[6])
                    a1 = float(params[7])
                    a2 = float(params[8])
                    a3 = float(params[9])
                    a4 = float(params[10])
                    if self.combined_bending_torsion_force is None:
                        self.combined_bending_torsion_force = mm.CustomCompoundBondForce(
                            4,
                            "k*sintheta0^3*sintheta1^3*(a0 + a1*cosphi + a2*cosphi^2 + a3*cosphi^3 + a4*cosphi^4); "
                            "sintheta0 = sin(angle(p1, p2, p3));"
                            "sintheta1 = sin(angle(p2, p3, p4));"
                            "cosphi = cos(dihedral(p1, p2, p3, p4));",
                        )
                        self.combined_bending_torsion_force.addPerBondParameter("k")
                        self.combined_bending_torsion_force.addPerBondParameter("a0")
                        self.combined_bending_torsion_force.addPerBondParameter("a1")
                        self.combined_bending_torsion_force.addPerBondParameter("a2")
                        self.combined_bending_torsion_force.addPerBondParameter("a3")
                        self.combined_bending_torsion_force.addPerBondParameter("a4")
                        sys.addForce(self.combined_bending_torsion_force)

                    self.combined_bending_torsion_force.addBond(
                        [
                            base_atom_index + atoms[0],
                            base_atom_index + atoms[1],
                            base_atom_index + atoms[2],
                            base_atom_index + atoms[3],
                        ],
                        [k, a0, a1, a2, a3, a4],
                    )
                else:
                    # RB Torsion
                    c = [float(x) for x in params[5:11]]
                    if any(x != 0 for x in c):
                        if self.rb_torsion_force is None:
                            self.rb_torsion_force = mm.RBTorsionForce()
                            sys.addForce(self.rb_torsion_force)
                        if dihedralType == "5":
                            # Convert Fourier coefficients to RB coefficients.
                            c = [
                                c[1] + 0.5 * (c[0] + c[2]),
                                0.5 * (-c[0] + 3 * c[2]),
                                -c[1] + 4 * c[3],
                                -2 * c[2],
                                -4 * c[3],
                                0,
                            ]
                        self.rb_torsion_force.addTorsion(
                            base_atom_index + atoms[0],
                            base_atom_index + atoms[1],
                            base_atom_index + atoms[2],
                            base_atom_index + atoms[3],
                            c[0],
                            c[1],
                            c[2],
                            c[3],
                            c[4],
                            c[5],
                        )

    def _add_cmap_to_system(self, sys, molecule_type, bonded_types, base_atom_index):
        cmap = None
        mapIndices = {}
        for fields in molecule_type.cmaps:
            atoms = [int(x) - 1 for x in fields[:5]]
            types = tuple(bonded_types[i] for i in atoms)
            if len(fields) >= 8 and len(fields) >= 8 + int(fields[6]) * int(fields[7]):
                params = fields
            elif types in self._cmapTypes:
                params = self._cmapTypes[types]
            elif types[::-1] in self._cmapTypes:
                params = self._cmapTypes[types[::-1]]
            else:
                raise ValueError(
                    "No parameters specified for cmap: "
                    + fields[0]
                    + ", "
                    + fields[1]
                    + ", "
                    + fields[2]
                    + ", "
                    + fields[3]
                    + ", "
                    + fields[4]
                )
            if cmap is None:
                cmap = mm.CMAPTorsionForce()
                sys.addForce(cmap)
            mapSize = int(params[6])
            if mapSize != int(params[7]):
                raise ValueError("Non-square CMAPs are not supported")
            map = []
            for i in range(mapSize):
                for j in range(mapSize):
                    map.append(
                        float(
                            params[
                                8
                                + mapSize * ((j + mapSize // 2) % mapSize)
                                + ((i + mapSize // 2) % mapSize)
                            ]
                        )
                    )
            map = tuple(map)
            if map not in mapIndices:
                mapIndices[map] = cmap.addMap(mapSize, map)
            cmap.addTorsion(
                mapIndices[map],
                base_atom_index + atoms[0],
                base_atom_index + atoms[1],
                base_atom_index + atoms[2],
                base_atom_index + atoms[3],
                base_atom_index + atoms[1],
                base_atom_index + atoms[2],
                base_atom_index + atoms[3],
                base_atom_index + atoms[4],
            )

    def _set_nonbonded_parameters(
        self,
        molecule_type,
        base_atom_index,
        atom_types,
        atom_type_map,
    ):
        charges = []
        for i, fields in enumerate(molecule_type.atoms):
            params = self._atom_types[fields[1]]
            if len(fields) > 6:
                q = float(fields[6])
            else:
                q = float(params[4])
            charges.append(q)
            index = base_atom_index + i
            atomType = atom_type_map[params[0]]
            self.nb_force.addParticle([atomType, q])
            # Add in the self term for the reaction field correction
            self.es_self_excl_force.addBond(index, index, [0.5 * q * q])

        pairExceptions = self._gen_pair_exceptions(
            molecule_type, base_atom_index, atom_types
        )
        exclusion_exceptions = self._gen_exclusion_exceptions(
            molecule_type, base_atom_index, atom_types
        )
        bondedExceptions = self._gen_bonded_exceptions(molecule_type, base_atom_index)
        constraint_exceptions = self._gen_constraint_exceptions(
            molecule_type, base_atom_index
        )
        # order matters
        exceptions = (
            pairExceptions
            + exclusion_exceptions
            + bondedExceptions
            + constraint_exceptions
        )
        return exceptions, charges

    def _gen_pair_exceptions(self, molecule_type, base_atom_index, atom_types):
        pairExceptions = []
        for fields in molecule_type.pairs:
            atoms = [int(x) - 1 for x in fields[:2]]
            types = tuple(atom_types[i] for i in atoms)
            if len(fields) >= 5:
                params = fields[3:5]
            elif types in self._pairTypes:
                params = self._pairTypes[types][3:5]
            elif types[::-1] in self._pairTypes:
                params = self._pairTypes[types[::-1]][3:5]
            else:
                # We'll use the automatically generated parameters
                continue

            type1_, q1 = self.nb_force.getParticleParameters(base_atom_index + atoms[0])
            type2_, q2 = self.nb_force.getParticleParameters(base_atom_index + atoms[1])
            pairExceptions.append(
                (
                    base_atom_index + atoms[0],
                    base_atom_index + atoms[1],
                    q1 * q2,
                    params[0],
                    params[1],
                )
            )
        return pairExceptions

    def _gen_exclusion_exceptions(self, molecule_type, base_atom_index, atom_types):
        exclusion_exceptions = []
        for fields in molecule_type.exclusions:
            atoms = [int(x) - 1 for x in fields]
            for atom in atoms[1:]:
                if atom > atoms[0]:
                    exclusion_exceptions.append(
                        (base_atom_index + atoms[0], base_atom_index + atom, 0, 0, 0)
                    )
        return exclusion_exceptions

    def _gen_bonded_exceptions(self, molecule_type, base_atom_index):
        bond_indices = []
        for fields in molecule_type.bonds:
            atoms = [int(x) - 1 for x in fields[:2]]
            bond_type = int(fields[2])
            # Type 6 bonds do not generate exclusions
            if bond_type == 1:
                bond_indices.append(
                    (base_atom_index + atoms[0], base_atom_index + atoms[1])
                )

        bondedExceptions = [(i, j, 0, 0, 0) for i, j in bond_indices]
        return bondedExceptions

    def _gen_constraint_exceptions(self, molecule_type, base_atom_index):
        constraint_indices = []
        for fields in molecule_type.constraints:
            atoms = [int(x) - 1 for x in fields[:2]]
            constraint_indices.append(
                (base_atom_index + atoms[0], base_atom_index + atoms[1])
            )

        constraint_exceptions = [(i, j, 0, 0, 0) for i, j in constraint_indices]
        return constraint_exceptions

    def _add_vsites_to_system(self, sys, molecule_type, base_atom_index):
        offset = base_atom_index - 1
        molecule_type.vsites.convert_com_to_linear(sys, offset)
        for index, site in molecule_type.vsites.iter():
            if isinstance(site, OutOfPlane):
                self._add_out_of_plane_vsite(sys, index, site, offset)
            elif isinstance(site, LinearSite):
                self._add_linear_vsite(sys, index, site, offset)
            elif isinstance(site, NormalizedInPlaneSite):
                self._add_normalized_in_plane_vsite(sys, index, site, offset)
            elif isinstance(site, NormalizedInPlaneTwoParticleSite):
                self._add_normalized_in_plane_two_particle_vsite(
                    sys, index, site, offset
                )
            else:
                raise RuntimeError(f"Unknown site type {type(site)}.")
            self._all_vsites.append(index + offset)

    def _add_normalized_in_plane_two_particle_vsite(self, sys, index, site, offset):
        vsite = mm.LocalCoordinatesSite(
            site.atom1 + offset,
            site.atom2 + offset,
            site.atom1 + offset,
            [1.0, 0.0, 0.0],
            [-0.5, 0.5, 0.0],
            [0.0, 0.0, 0.0],
            [site.a, 0.0, 0.0],
        )
        sys.setVirtualSite(index + offset, vsite)

    def _add_normalized_in_plane_vsite(self, sys, index, site, offset):
        vsite = mm.LocalCoordinatesSite(
            site.atom1 + offset,
            site.atom2 + offset,
            site.atom3 + offset,
            [1.0, 0.0, 0.0],
            [-1.0, 1.0 - site.a, site.a],
            [0.0, 0.0, 0.0],
            [site.d, 0.0, 0.0],
        )
        sys.setVirtualSite(index + offset, vsite)

    def _add_out_of_plane_vsite(self, sys, index, site, offset):
        vsite = mm.OutOfPlaneSite(
            site.atom1 + offset,
            site.atom2 + offset,
            site.atom3 + offset,
            site.a,
            site.b,
            site.c,
        )
        sys.setVirtualSite(index + offset, vsite)

    def _add_linear_vsite(self, sys, index, site, offset):
        n = len(site.atom_weights)
        if n == 1:
            self._add_one_particle_vsite(sys, index, site, offset)
        elif n == 2:
            self._add_two_particle_vsite(sys, index, site, offset)
        elif n == 3:
            self._add_three_particle_vsite(sys, index, site, offset)
        else:
            self._add_n_particle_vsite(sys, index, site, offset)

    def _add_one_particle_vsite(self, sys, index, site, offset):
        atoms = []
        weights = []
        for atom, weight in site.atom_weights.items():
            atoms.append(atom + offset)
            weights.append(weight)
        assert len(atoms) == 1

        # We need at least two atoms, so we'll use atom 0 as
        # dummy with a weight of 0.
        atoms.append(0)
        weights.append(0.0)

        # These are dummy weights that specify a local coordinate system
        # that openmm supports, but we don't use.
        x_weights = [0.0, 0.0]
        y_weights = [0.0, 0.0]

        vsite = mm.LocalCoordinatesSite(
            atoms, weights, x_weights, y_weights, [0.0, 0.0, 0.0]
        )
        sys.setVirtualSite(index + offset, vsite)

        # Add exclusions
        self.nb_force.addExclusion(index + offset, atoms[0])

    def _add_two_particle_vsite(self, sys, index, site, offset):
        atoms = []
        weights = []
        for atom, weight in site.atom_weights.items():
            atoms.append(atom + offset)
            weights.append(weight)
        assert len(atoms) == 2

        vsite = mm.TwoParticleAverageSite(atoms[0], atoms[1], weights[0], weights[1])
        sys.setVirtualSite(index + offset, vsite)

    def _add_three_particle_vsite(self, sys, index, site, offset):
        atoms = []
        weights = []
        for atom, weight in site.atom_weights.items():
            atoms.append(atom + offset)
            weights.append(weight)
        assert len(atoms) == 3

        vsite = mm.ThreeParticleAverageSite(
            atoms[0],
            atoms[1],
            atoms[2],
            weights[0],
            weights[1],
            weights[2],
        )
        sys.setVirtualSite(index + offset, vsite)

    def _add_n_particle_vsite(self, sys, index, site, offset):
        atoms = []
        weights = []
        for atom, weight in site.atom_weights.items():
            atoms.append(atom + offset)
            weights.append(weight)
        n_atoms = len(atoms)

        # These are dummy weights that specify a local coordinate system
        # that openmm supports, but we don't use.
        x_weights = [0.0] * n_atoms
        y_weights = [0.0] * n_atoms
        vsite = mm.LocalCoordinatesSite(
            atoms, weights, x_weights, y_weights, [0.0, 0.0, 0.0]
        )
        sys.setVirtualSite(index + offset, vsite)

    def create_system(
        self,
        nonbonded_cutoff=1.1 * unit.nanometer,
        remove_com_motion=True,
    ):
        """Construct an OpenMM System representing the topology described by this
        top file.

        Parameters
        ----------
        nonbonded_cutoff : distance=1*nanometer
            The cutoff distance to use for nonbonded interactions
        remove_com_motion : boolean=True
            If true, a CMMotionRemover will be added to the System
        Returns
        -------
        System
             the newly created System
        """
        sys = mm.System()
        box_vectors = self.topology.getPeriodicBoxVectors()
        if box_vectors is not None:
            sys.setDefaultPeriodicBoxVectors(*box_vectors)
        else:
            raise ValueError("periodicBoxVectors must be set")

        # Custom force to handle martini nonbonded terms
        # nbfix-like terms
        self.nb_force = mm.CustomNonbondedForce(
            "step(rcut-r)*(LJ - corr + ES);"
            "LJ = (C12(type1, type2) / r^12 - C6(type1, type2) / r^6);"
            "corr = (C12(type1, type2) / rcut^12 - C6(type1, type2) / rcut^6);"
            "ES = f/epsilon_r*q1*q2 * (1/r + krf * r^2 - crf);"
            "crf = 1 / rcut + krf * rcut^2;"
            "krf = 1 / (2 * rcut^3);"
            f"epsilon_r = {self.epsilon_r};"
            "f = 138.935458;"
            f"rcut={nonbonded_cutoff.value_in_unit(unit.nanometers)};"
        )
        self.nb_force.addPerParticleParameter("type")
        self.nb_force.addPerParticleParameter("q")
        self.nb_force.setNonbondedMethod(mm.CustomNonbondedForce.CutoffPeriodic)
        self.nb_force.setCutoffDistance(nonbonded_cutoff.value_in_unit(unit.nanometer))
        sys.addForce(self.nb_force)

        # custom non-bonded force to add in the electrostatic
        # terms for self and excluded interactions
        self.es_self_excl_force = mm.CustomBondForce(
            f"step(rcut-r) * ES;"
            f"ES = f*qprod/epsilon_r * (krf * r^2 - crf);"
            f"crf = 1 / rcut + krf * rcut^2;"
            f"krf = 1 / (2 * rcut^3);"
            f"epsilon_r = {self.epsilon_r};"
            f"f = 138.935458;"
            f"rcut={nonbonded_cutoff.value_in_unit(unit.nanometers)};"
        )
        self.es_self_excl_force.addPerBondParameter("qprod")
        sys.addForce(self.es_self_excl_force)

        # custom non-bonded force for the ES exceptions
        self.es_except_force = mm.CustomBondForce(
            f"step(rcut-r) * ES;"
            f"ES = f*qprod/epsilon_r * (1/r + krf * r^2 - crf);"
            f"crf = 1 / rcut + krf * rcut^2;"
            f"krf = 1 / (2 * rcut^3);"
            f"epsilon_r = {self.epsilon_r};"
            f"f = 138.935458;"
            f"rcut={nonbonded_cutoff.value_in_unit(unit.nanometers)};"
        )
        self.es_except_force.addPerBondParameter("qprod")
        sys.addForce(self.es_except_force)

        # custom bonded force to handle exceptions
        self.lj_except_force = mm.CustomBondForce(
            "step(rcut-r) * (energy - corr);"
            "energy = (C12/r^12 - C6/r^6);"
            "corr = (C12/rcut^12 - C6/rcut^6);"
            f"rcut={nonbonded_cutoff.value_in_unit(unit.nanometers)};"
        )
        self.lj_except_force.addPerBondParameter("C12")
        self.lj_except_force.addPerBondParameter("C6")
        sys.addForce(self.lj_except_force)

        all_exceptions = []
        all_charges = []

        (
            dihedral_type_table,
            wildcard_dihedral_types,
        ) = self._build_dihedral_lookup_table()

        # build a lookup table mapping atom types into integer indices
        # that will later be used to lookup LJ combinations
        used_atom_types = set()
        for molecule_name, _ in self._molecules:
            molecule_type = self._moleculeTypes[molecule_name]
            for atom in molecule_type.atoms:
                used_atom_types.add(atom[1])
        atom_type_map = {k: i for i, k in enumerate(sorted(used_atom_types))}

        # now we need to setup our tables of C6 and C12
        n_types = len(atom_type_map)
        C6 = []
        C12 = []
        for type_i in atom_type_map:
            for type_j in atom_type_map:
                type_i_sorted, type_j_sorted = sorted([type_i, type_j])
                if (type_i_sorted, type_j_sorted) in self._nonbond_types:
                    # Parameters were specified for this pair of types,
                    # so use them.
                    params = self._nonbond_types[(type_i_sorted, type_j_sorted)]
                    if self._use_sigma_eps:
                        sigma = float(params[3])
                        eps = float(params[4])
                        c6 = 4 * eps * sigma ** 6
                        c12 = 4 * eps * sigma ** 12
                    else:
                        c6 = float(params[3])
                        c12 = float(params[4])
                else:
                    # Parameters were not specified for this pair type,
                    # so calculate using combination rules.
                    params_i = self._atom_types[type_i]
                    v_i = float(params_i[6])
                    w_i = float(params_i[7])
                    params_j = self._atom_types[type_j]
                    v_j = float(params_j[6])
                    w_j = float(params_j[7])

                    if self._use_sigma_eps:
                        sigma = 0.5 * (v_i + v_j)
                        eps = math.sqrt(w_i * w_j)
                        c6 = 4 * eps * sigma ** 6
                        c12 = 4 * eps * sigma ** 12
                    else:
                        c6 = math.sqrt(v_i * v_j)
                        c12 = math.sqrt(w_i * w_j)

                C6.append(c6)
                C12.append(c12)
        self.nb_force.addTabulatedFunction(
            "C6", mm.Discrete2DFunction(n_types, n_types, C6)
        )
        self.nb_force.addTabulatedFunction(
            "C12", mm.Discrete2DFunction(n_types, n_types, C12)
        )

        # Loop over molecules and create the specified number of each type.
        for molecule_name, molecule_count in self._molecules:
            molecule_type = self._moleculeTypes[molecule_name]
            for i in range(molecule_count):
                base_atom_index = sys.getNumParticles()
                atom_types = [atom[1] for atom in molecule_type.atoms]
                try:
                    bonded_types = [self._atom_types[t][1] for t in atom_types]
                except KeyError as e:
                    raise ValueError("Unknown atom type: " + e.message)
                bonded_types = [
                    b if b is not None else a for a, b in zip(atom_types, bonded_types)
                ]

                # add bonded parameters
                self._add_atoms_to_system(sys, molecule_type)
                self._add_bonds_to_system(
                    sys, molecule_type, bonded_types, base_atom_index
                )
                self._add_angles_to_system(
                    sys, molecule_type, bonded_types, base_atom_index
                )
                self._add_torsions_to_system(
                    sys,
                    molecule_type,
                    bonded_types,
                    dihedral_type_table,
                    wildcard_dihedral_types,
                    base_atom_index,
                )
                self._add_cmap_to_system(
                    sys, molecule_type, bonded_types, base_atom_index
                )
                self._add_vsites_to_system(sys, molecule_type, base_atom_index)

                exceptions, charges = self._set_nonbonded_parameters(
                    molecule_type,
                    base_atom_index,
                    atom_types,
                    atom_type_map,
                )
                all_exceptions.extend(exceptions)
                all_charges.extend(charges)

                # Add explicitly specified constraints.
                for fields in molecule_type.constraints:
                    atoms = [int(x) - 1 for x in fields[:2]]
                    constraint_type = int(fields[2])
                    assert constraint_type == 1
                    length = float(fields[3])
                    sys.addConstraint(
                        base_atom_index + atoms[0], base_atom_index + atoms[1], length
                    )

        if all_exceptions:
            # build a map of unique exceptions
            # process in order, so that later entries trump earlier ones
            except_map = defaultdict(list)
            for exception in all_exceptions:
                i, j, q, c6, c12 = exception
                if i < j:
                    except_map[(i, j)] = [q, c6, c12]
                else:
                    except_map[(j, i)] = [q, c6, c12]

            # add in all of the exclusions
            for i, j in except_map:
                # Remove i,j from nonbonded interactions for all exceptions / exclusions
                self.nb_force.addExclusion(i, j)
            
                # Handle electrostatic exceptions / exclusions.
                # We're going to assume that q==0 means that this was an
                # exclusion.
                # In this for loop, we have not give the q, c6, c12 correspond to i and j.
                q = float(except_map[(i, j)][0])
                c6 = float(except_map[(i, j)][1])
                c12 = float(except_map[(i, j)][0])
                if q == 0:
                    # In this case, we still need to add in the reaction field correction
                    # term.
                    qprod = all_charges[i] * all_charges[j]
                    # Don't bother adding the interaction if either particle has zero charge.
                    if qprod != 0.0:
                        self.es_self_excl_force.addBond(i, j, [qprod])
                # If q !=0, then this is an exception.
                else:
                    # We don't bother adding the interaction if either the charge product is zero
                    self.es_except_force.addBond(i, j, [q])
                # Now we'll add the LJ exceptions. We don't bother
                # adding the interaction if the combined LJ parameters are zero.
                if c6 != 0 and c12 != 0:
                # As in lj_except_force we first add the parameter C12 and then C6
                    self.lj_except_force.addBond(i, j, [c12, c6])

        if remove_com_motion:
            sys.addForce(mm.CMMotionRemover())
        return sys


def _get_default_gromacs_include_dir():
    """Find the location where gromacs #include files are referenced from, by
    searching for (1) gromacs environment variables, (2) for the gromacs binary
    'pdb2gmx' or 'gmx' in the PATH, or (3) just using the default gromacs
    install location, /usr/local/gromacs/share/gromacs/top"""
    if "GMXDATA" in os.environ:
        return os.path.join(os.environ["GMXDATA"], "top")
    if "GMXBIN" in os.environ:
        return os.path.abspath(
            os.path.join(os.environ["GMXBIN"], "..", "share", "gromacs", "top")
        )

    pdb2gmx_path = distutils.spawn.find_executable("pdb2gmx")
    if pdb2gmx_path is not None:
        return os.path.abspath(
            os.path.join(os.path.dirname(pdb2gmx_path), "..", "share", "gromacs", "top")
        )
    else:
        gmx_path = distutils.spawn.find_executable("gmx")
        if gmx_path is not None:
            return os.path.abspath(
                os.path.join(os.path.dirname(gmx_path), "..", "share", "gromacs", "top")
            )

    return "/usr/local/gromacs/share/gromacs/top"
