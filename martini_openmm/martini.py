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
from __future__ import absolute_import

__author__ = "Justin MacCallum"
__version__ = "0.1"

from simtk.openmm.app import Topology
from simtk.openmm.app import PDBFile
from simtk.openmm.app import forcefield as ff
from simtk.openmm.app import element as elem
from simtk.openmm.app import amberprmtopfile as prmtop
import simtk.unit as unit
import simtk.openmm as mm
import math
import os
import re
import distutils.spawn
from collections import OrderedDict, defaultdict
from itertools import combinations, combinations_with_replacement

HBonds = ff.HBonds
AllBonds = ff.AllBonds
HAngles = ff.HAngles

OBC2 = prmtop.OBC2

novarcharre = re.compile(r"\W")


def _find_all_instances_in_string(string, substr):
    """ Find indices of all instances of substr in string """
    indices = []
    idx = string.find(substr, 0)
    while idx > -1:
        indices.append(idx)
        idx = string.find(substr, idx + 1)
    return indices


def _replace_defines(line, defines):
    """ Replaces defined tokens in a given line """
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


class GromacsMartiniV2TopFile(object):
    """GromacsMartiniV2TopFile parses a Gromacs top file and constructs a Topology and (optionally) an OpenMM System from it."""

    class _MoleculeType(object):
        """Inner class to store information about a molecule type."""

        def __init__(self):
            self.atoms = []
            self.bonds = []
            self.angles = []
            self.dihedrals = []
            self.exclusions = []
            self.pairs = []
            self.constraints = []
            self.cmaps = []
            self.vsites2 = []
            self.has_virtual_sites = False
            self.has_nbfix_terms = False

    def _processFile(self, file):
        append = ""
        for line in open(file):
            if line.strip().endswith("\\"):
                append = "%s %s" % (append, line[: line.rfind("\\")])
            else:
                self._processLine(append + " " + line, file)
                append = ""

    def _processLine(self, line, file):
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
                        self._processFile(file)
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
                del (self._ifStack[-1])
                del (self._elseStack[-1])
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
                self._processDefaults(line)
            elif self._currentCategory == "moleculetype":
                self._processMoleculeType(line)
            elif self._currentCategory == "molecules":
                self._processMolecule(line)
            elif self._currentCategory == "atoms":
                self._processAtom(line)
            elif self._currentCategory == "bonds":
                self._processBond(line)
            elif self._currentCategory == "angles":
                self._processAngle(line)
            elif self._currentCategory == "dihedrals":
                self._processDihedral(line)
            elif self._currentCategory == "exclusions":
                self._processExclusion(line)
            elif self._currentCategory == "pairs":
                self._processPair(line)
            elif self._currentCategory == "constraints":
                self._processConstraint(line)
            elif self._currentCategory == "cmap":
                self._processCmap(line)
            elif self._currentCategory == "atomtypes":
                self._processAtomType(line)
            elif self._currentCategory == "bondtypes":
                self._processBondType(line)
            elif self._currentCategory == "angletypes":
                self._processAngleType(line)
            elif self._currentCategory == "dihedraltypes":
                self._processDihedralType(line)
            elif self._currentCategory == "implicit_genborn_params":
                self._processImplicitType(line)
            elif self._currentCategory == "pairtypes":
                self._processPairType(line)
            elif self._currentCategory == "cmaptypes":
                self._processCmapType(line)
            elif self._currentCategory == "nonbond_params":
                self._processNonbondType(line)
            elif self._currentCategory == "virtual_sites2":
                self._processVirtualSites2(line)
            elif self._currentCategory.startswith("virtual_sites"):
                if self._currentMoleculeType is None:
                    raise ValueError(
                        "Found %s before [ moleculetype ]" % self._currentCategory
                    )
                self._currentMoleculeType.has_virtual_sites = True

    def _processDefaults(self, line):
        """Process the [ defaults ] line."""

        fields = line.split()
        if len(fields) > 2:
            raise ValueError("Too many fields in [ defaults ] line: " + line)
        if fields[0] != "1":
            raise ValueError("Unsupported nonbonded type: " + fields[0])
        if fields[1] != "1":
            raise ValueError("Unsupported combination rule: " + fields[1])
        self._defaults = fields

    def _processMoleculeType(self, line):
        """Process a line in the [ moleculetypes ] category."""
        fields = line.split()
        if len(fields) < 1:
            raise ValueError("Too few fields in [ moleculetypes ] line: " + line)
        type = GromacsMartiniV2TopFile._MoleculeType()
        self._moleculeTypes[fields[0]] = type
        self._currentMoleculeType = type

    def _processMolecule(self, line):
        """Process a line in the [ molecules ] category."""
        fields = line.split()
        if len(fields) < 2:
            raise ValueError("Too few fields in [ molecules ] line: " + line)
        self._molecules.append((fields[0], int(fields[1])))

    def _processAtom(self, line):
        """Process a line in the [ atoms ] category."""
        if self._currentMoleculeType is None:
            raise ValueError("Found [ atoms ] section before [ moleculetype ]")
        fields = line.split()
        if len(fields) < 5:
            raise ValueError("Too few fields in [ atoms ] line: " + line)
        self._currentMoleculeType.atoms.append(fields)

    def _processBond(self, line):
        """Process a line in the [ bonds ] category."""
        if self._currentMoleculeType is None:
            raise ValueError("Found [ bonds ] section before [ moleculetype ]")
        fields = line.split()
        if len(fields) < 3:
            raise ValueError("Too few fields in [ bonds ] line: " + line)
        if fields[2] != "1":
            raise ValueError("Unsupported function type in [ bonds ] line: " + line)
        self._currentMoleculeType.bonds.append(fields)

    def _processAngle(self, line):
        """Process a line in the [ angles ] category."""
        if self._currentMoleculeType is None:
            raise ValueError("Found [ angles ] section before [ moleculetype ]")
        fields = line.split()
        if len(fields) < 4:
            raise ValueError("Too few fields in [ angles ] line: " + line)
        if fields[3] not in ("2"):
            raise ValueError(
                "Unsupported (nonG96) function type in [ angles ] line: " + line
            )
        self._currentMoleculeType.angles.append(fields)

    def _processDihedral(self, line):
        """Process a line in the [ dihedrals ] category."""
        if self._currentMoleculeType is None:
            raise ValueError("Found [ dihedrals ] section before [ moleculetype ]")
        fields = line.split()
        if len(fields) < 5:
            raise ValueError("Too few fields in [ dihedrals ] line: " + line)
        if fields[4] not in ("1", "2", "3", "4", "5", "9"):
            raise ValueError("Unsupported function type in [ dihedrals ] line: " + line)
        self._currentMoleculeType.dihedrals.append(fields)

    def _processExclusion(self, line):
        """Process a line in the [ exclusions ] category."""
        if self._currentMoleculeType is None:
            raise ValueError("Found [ exclusions ] section before [ moleculetype ]")
        fields = line.split()
        if len(fields) < 2:
            raise ValueError("Too few fields in [ exclusions ] line: " + line)
        self._currentMoleculeType.exclusions.append(fields)

    def _processPair(self, line):
        """Process a line in the [ pairs ] category."""
        if self._currentMoleculeType is None:
            raise ValueError("Found [ pairs ] section before [ moleculetype ]")
        fields = line.split()
        if len(fields) < 3:
            raise ValueError("Too few fields in [ pairs ] line: " + line)
        if fields[2] != "1":
            raise ValueError("Unsupported function type in [ pairs ] line: " + line)
        self._currentMoleculeType.pairs.append(fields)

    def _processConstraint(self, line):
        """Process a line in the [ constraints ] category."""
        if self._currentMoleculeType is None:
            raise ValueError("Found [ constraints ] section before [ moleculetype ]")
        fields = line.split()
        if len(fields) < 4:
            raise ValueError("Too few fields in [ constraints ] line: " + line)
        self._currentMoleculeType.constraints.append(fields)

    def _processCmap(self, line):
        """Process a line in the [ cmaps ] category."""
        if self._currentMoleculeType is None:
            raise ValueError("Found [ cmap ] section before [ moleculetype ]")
        fields = line.split()
        if len(fields) < 6:
            raise ValueError("Too few fields in [ cmap ] line: " + line)
        self._currentMoleculeType.cmaps.append(fields)

    def _processAtomType(self, line):
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
        self._atomTypes[fields[0]] = fields

    def _processBondType(self, line):
        """Process a line in the [ bondtypes ] category."""
        fields = line.split()
        if len(fields) < 5:
            raise ValueError("Too few fields in [ bondtypes ] line: " + line)
        if fields[2] != "1":
            raise ValueError("Unsupported function type in [ bondtypes ] line: " + line)
        self._bondTypes[tuple(fields[:2])] = fields

    def _processAngleType(self, line):
        """Process a line in the [ angletypes ] category."""
        fields = line.split()
        if len(fields) < 6:
            raise ValueError("Too few fields in [ angletypes ] line: " + line)
        if fields[3] not in ("1", "5"):
            raise ValueError(
                "Unsupported function type in [ angletypes ] line: " + line
            )
        self._angleTypes[tuple(fields[:3])] = fields

    def _processDihedralType(self, line):
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

    def _processImplicitType(self, line):
        """Process a line in the [ implicit_genborn_params ] category."""
        fields = line.split()
        if len(fields) < 6:
            raise ValueError(
                "Too few fields in [ implicit_genborn_params ] line: " + line
            )
        self._implicitTypes[fields[0]] = fields

    def _processPairType(self, line):
        """Process a line in the [ pairtypes ] category."""
        fields = line.split()
        if len(fields) < 5:
            raise ValueError("Too few fields in [ pairtypes] line: " + line)
        if fields[2] != "1":
            raise ValueError("Unsupported function type in [ pairtypes ] line: " + line)
        self._pairTypes[tuple(fields[:2])] = fields

    def _processCmapType(self, line):
        """Process a line in the [ cmaptypes ] category."""
        fields = line.split()
        if len(fields) < 8 or len(fields) < 8 + int(fields[6]) * int(fields[7]):
            raise ValueError("Too few fields in [ cmaptypes ] line: " + line)
        if fields[5] != "1":
            raise ValueError("Unsupported function type in [ cmaptypes ] line: " + line)
        self._cmapTypes[tuple(fields[:5])] = fields

    def _processNonbondType(self, line):
        """Process a line in the [ nonbond_params ] category."""
        fields = line.split()
        if len(fields) < 5:
            raise ValueError("Too few fields in [ nonbond_params ] line: " + line)
        if fields[2] != "1":
            raise ValueError(
                "Unsupported function type in [ nonbond_params ] line: " + line
            )
        self._nonbondTypes[tuple(sorted(fields[:2]))] = fields

    def _processVirtualSites2(self, line):
        """Process a line in the [ virtual_sites2 ] category."""
        fields = line.split()
        if len(fields) < 5:
            raise ValueError("Too few fields in [ virtual_sites2 ] line: " + line)
        self._currentMoleculeType.vsites2.append(fields[:5])

    def _buildDihLookupTable(self):
        dihedralTypeTable = {}
        for key in self._dihedralTypes:
            if key[1] != "X" and key[2] != "X":
                if (key[1], key[2]) not in dihedralTypeTable:
                    dihedralTypeTable[(key[1], key[2])] = []
                dihedralTypeTable[(key[1], key[2])].append(key)
                if (key[2], key[1]) not in dihedralTypeTable:
                    dihedralTypeTable[(key[2], key[1])] = []
                dihedralTypeTable[(key[2], key[1])].append(key)
        wildcardDihedralTypes = []
        for key in self._dihedralTypes:
            if key[1] == "X" or key[2] == "X":
                wildcardDihedralTypes.append(key)
                for types in dihedralTypeTable.values():
                    types.append(key)
        return dihedralTypeTable, wildcardDihedralTypes

    def _addAtomsToSystem(self, sys, moleculeType):
        for fields in moleculeType.atoms:
            if len(fields) >= 8:
                mass = float(fields[7])
            else:
                mass = float(self._atomTypes[fields[1]][3])
            sys.addParticle(mass)

    def _addBondsToSystem(self, sys, moleculeType, bondedTypes, baseAtomIndex):
        bonds = None
        for fields in moleculeType.bonds:
            atoms = [int(x) - 1 for x in fields[:2]]
            types = tuple(bondedTypes[i] for i in atoms)
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

            if bonds is None:
                bonds = mm.HarmonicBondForce()
                sys.addForce(bonds)
            bonds.addBond(
                baseAtomIndex + atoms[0],
                baseAtomIndex + atoms[1],
                length,
                float(params[1]),
            )

    def _addAnglesToSystem(self, sys, moleculeType, bondedTypes, baseAtomIndex):
        angles_harmonic = None
        angles_gmx = None
        degToRad = math.pi / 180

        for fields in moleculeType.angles:
            atoms = [int(x) - 1 for x in fields[:3]]
            types = tuple(bondedTypes[i] for i in atoms)
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

            # TODO This part is broken
            # Need to implment gmx angle types
            theta = float(params[0]) * degToRad
            if int(fields[3]) == 2:
                if angles_gmx is None:
                    angles_gmx = mm.CustomAngleForce("0.5*k*(cos(theta)-cos(theta0))^2")
                    angles_gmx.addPerAngleParameter("theta0")
                    angles_gmx.addPerAngleParameter("k")
                    sys.addForce(angles_gmx)
                angles_gmx.addAngle(
                    baseAtomIndex + atoms[0],
                    baseAtomIndex + atoms[1],
                    baseAtomIndex + atoms[2],
                    [theta, float(params[1])],
                )

            elif int(fields[3]) == 1:
                if angles_harmonic is None:
                    angles_harmonic = mm.HarmonicAngleForce()
                    sys.addForce(angles_harmonic)
                angles_harmonic.addAngle(
                    baseAtomIndex + atoms[0],
                    baseAtomIndex + atoms[1],
                    baseAtomIndex + atoms[2],
                    theta,
                    float(params[1]),
                )

            else:
                raise RuntimeError()

    def _addTorsionToSystem(
        self,
        sys,
        moleculeType,
        bondedTypes,
        dihedralTypeTable,
        wildcardDihedralTypes,
        baseAtomIndex,
    ):
        periodic = None
        harmonicTorsion = None
        rb = None
        degToRad = math.pi / 180

        for fields in moleculeType.dihedrals:
            atoms = [int(x) - 1 for x in fields[:4]]
            types = tuple(bondedTypes[i] for i in atoms)
            dihedralType = fields[4]
            reversedTypes = types[::-1] + (dihedralType,)
            types = types + (dihedralType,)
            if (
                (dihedralType in ("1", "4", "5", "9") and len(fields) > 7)
                or (dihedralType == "3" and len(fields) > 10)
                or (dihedralType == "2" and len(fields) > 6)
            ):
                paramsList = [fields]
            else:
                # Look for a matching dihedral type.
                paramsList = None
                if (types[1], types[2]) in dihedralTypeTable:
                    dihedralTypes = dihedralTypeTable[(types[1], types[2])]
                else:
                    dihedralTypes = wildcardDihedralTypes
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
                        if periodic is None:
                            periodic = mm.PeriodicTorsionForce()
                            sys.addForce(periodic)
                        periodic.addTorsion(
                            baseAtomIndex + atoms[0],
                            baseAtomIndex + atoms[1],
                            baseAtomIndex + atoms[2],
                            baseAtomIndex + atoms[3],
                            int(float(params[7])),
                            float(params[5]) * degToRad,
                            k,
                        )
                elif dihedralType == "2":
                    # Harmonic torsion
                    k = float(params[6])
                    phi0 = float(params[5])
                    if k != 0:
                        if harmonicTorsion is None:
                            harmonicTorsion = mm.CustomTorsionForce(
                                "0.5*k*(thetap-theta0)^2; thetap = step(-(theta-theta0+pi))*2*pi+theta+step(theta-theta0-pi)*(-2*pi); pi = %.15g"
                                % math.pi
                            )
                            harmonicTorsion.addPerTorsionParameter("theta0")
                            harmonicTorsion.addPerTorsionParameter("k")
                            sys.addForce(harmonicTorsion)
                        # map phi0 into correct space
                        phi0 = phi0 - 360 if phi0 > 180 else phi0
                        harmonicTorsion.addTorsion(
                            baseAtomIndex + atoms[0],
                            baseAtomIndex + atoms[1],
                            baseAtomIndex + atoms[2],
                            baseAtomIndex + atoms[3],
                            (phi0 * degToRad, k),
                        )
                else:
                    # RB Torsion
                    c = [float(x) for x in params[5:11]]
                    if any(x != 0 for x in c):
                        if rb is None:
                            rb = mm.RBTorsionForce()
                            sys.addForce(rb)
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
                        rb.addTorsion(
                            baseAtomIndex + atoms[0],
                            baseAtomIndex + atoms[1],
                            baseAtomIndex + atoms[2],
                            baseAtomIndex + atoms[3],
                            c[0],
                            c[1],
                            c[2],
                            c[3],
                            c[4],
                            c[5],
                        )

    def _addCmapToSystem(self, sys, moleculeType, bondedTypes, baseAtomIndex):
        cmap = None
        mapIndices = {}
        for fields in moleculeType.cmaps:
            atoms = [int(x) - 1 for x in fields[:5]]
            types = tuple(bondedTypes[i] for i in atoms)
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
                baseAtomIndex + atoms[0],
                baseAtomIndex + atoms[1],
                baseAtomIndex + atoms[2],
                baseAtomIndex + atoms[3],
                baseAtomIndex + atoms[1],
                baseAtomIndex + atoms[2],
                baseAtomIndex + atoms[3],
                baseAtomIndex + atoms[4],
            )

    def _setnonbondedParams(
        self, nb, elecSelfForce, lj, moleculeType, baseAtomIndex, atomTypes, atomTypeMap
    ):
        charges = []
        for i, fields in enumerate(moleculeType.atoms):
            params = self._atomTypes[fields[1]]
            if len(fields) > 6:
                q = float(fields[6])
            else:
                q = float(params[4])
            charges.append(q)
            nb.addParticle([q])
            index = baseAtomIndex + i
            elecSelfForce.addBond(index, index, [0.5 * q * q])
            atomType = atomTypeMap[params[0]]
            lj.addParticle([atomType])

        pairExceptions = self._genPairExceptions(
            nb, moleculeType, baseAtomIndex, atomTypes
        )
        exclusionExceptions = self._genExclusionExceptions(
            nb, moleculeType, baseAtomIndex, atomTypes
        )
        bondedExceptions = self._genBondedExceptions(moleculeType, baseAtomIndex)
        constraintExceptions = self._genConstraintExceptions(
            moleculeType, baseAtomIndex
        )
        # order matters
        exceptions = (
            pairExceptions
            + exclusionExceptions
            + bondedExceptions
            + constraintExceptions
        )
        return exceptions, charges

    def _genPairExceptions(self, nb, moleculeType, baseAtomIndex, atomTypes):
        pairExceptions = []
        for fields in moleculeType.pairs:
            atoms = [int(x) - 1 for x in fields[:2]]
            types = tuple(atomTypes[i] for i in atoms)
            if len(fields) >= 5:
                params = fields[3:5]
            elif types in self._pairTypes:
                params = self._pairTypes[types][3:5]
            elif types[::-1] in self._pairTypes:
                params = self._pairTypes[types[::-1]][3:5]
            else:
                # We'll use the automatically generated parameters
                continue

            atom1params = nb.getParticleParameters(baseAtomIndex + atoms[0])
            atom2params = nb.getParticleParameters(baseAtomIndex + atoms[1])
            pairExceptions.append(
                (
                    baseAtomIndex + atoms[0],
                    baseAtomIndex + atoms[1],
                    atom1params[0] * atom2params[0],
                    params[0],
                    params[1],
                )
            )
        return pairExceptions

    def _genExclusionExceptions(self, nb, moleculeType, baseAtomIndex, atomTypes):
        exclusionExceptions = []
        for fields in moleculeType.exclusions:
            atoms = [int(x) - 1 for x in fields]
            for atom in atoms[1:]:
                if atom > atoms[0]:
                    exclusionExceptions.append(
                        (baseAtomIndex + atoms[0], baseAtomIndex + atom, 0, 0, 0)
                    )
        return exclusionExceptions

    def _genBondedExceptions(self, moleculeType, baseAtomIndex):
        bondIndices = []
        for fields in moleculeType.bonds:
            atoms = [int(x) - 1 for x in fields[:2]]
            bondIndices.append((baseAtomIndex + atoms[0], baseAtomIndex + atoms[1]))

        bondedExceptions = [(i, j, 0, 0, 0) for i, j in bondIndices]
        return bondedExceptions

    def _genConstraintExceptions(self, moleculeType, baseAtomIndex):
        constraintIndices = []
        for fields in moleculeType.constraints:
            atoms = [int(x) - 1 for x in fields[:2]]
            constraintIndices.append(
                (baseAtomIndex + atoms[0], baseAtomIndex + atoms[1])
            )

        constraintExceptions = [(i, j, 0, 0, 0) for i, j in constraintIndices]
        return constraintExceptions

    def __init__(
        self,
        file,
        periodicBoxVectors=None,
        unitCellDimensions=None,
        includeDir=None,
        defines=None,
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
        if includeDir is None:
            includeDir = _defaultGromacsIncludeDir()
        self._includeDirs = (os.path.dirname(file), includeDir)
        # Most of the gromacs water itp files for different forcefields,
        # unless the preprocessor #define FLEXIBLE is given, don't define
        # bonds between the water hydrogen and oxygens, but only give the
        # constraint distances and exclusions.
        self._defines = OrderedDict()
        self._defines["FLEXIBLE"] = True
        self._genpairs = True
        if defines is not None:
            for define, value in defines.iteritems():
                self._defines[define] = value

        # Parse the file.

        self._currentCategory = None
        self._ifStack = []
        self._elseStack = []
        self._moleculeTypes = {}
        self._molecules = []
        self._currentMoleculeType = None
        self._atomTypes = {}
        self._bondTypes = {}
        self._angleTypes = {}
        self._dihedralTypes = {}
        self._implicitTypes = {}
        self._pairTypes = {}
        self._cmapTypes = {}
        self._nonbondTypes = {}
        self._processFile(file)

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
        for moleculeName, moleculeCount in self._molecules:
            if moleculeName not in self._moleculeTypes:
                raise ValueError("Unknown molecule type: " + moleculeName)
            moleculeType = self._moleculeTypes[moleculeName]
            if moleculeCount > 0 and moleculeType.has_virtual_sites:
                raise ValueError("Virtual sites not yet supported by Gromacs parsers")

            # Create the specified number of molecules of this type.
            for i in range(moleculeCount):
                atoms = []
                lastResidue = None
                c = top.addChain()
                for index, fields in enumerate(moleculeType.atoms):
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

                for fields in moleculeType.bonds:
                    top.addBond(atoms[int(fields[0]) - 1], atoms[int(fields[1]) - 1])

    def createSystem(
        self,
        nonbondedMethod=ff.NoCutoff,
        nonbondedCutoff=1.1 * unit.nanometer,
        epsilon_r=15.0,
        removeCMMotion=True,
    ):
        """Construct an OpenMM System representing the topology described by this
        top file.

        Parameters
        ----------
        nonbondedMethod : object=NoCutoff
            The method to use for nonbonded interactions.  Allowed values are
            NoCutoff, CutoffNonPeriodic, CutoffPeriodic, Ewald, PME, or LJPME.
        nonbondedCutoff : distance=1*nanometer
            The cutoff distance to use for nonbonded interactions
        removeCMMotion : boolean=True
            If true, a CMMotionRemover will be added to the System
        Returns
        -------
        System
             the newly created System
        """
        sys = mm.System()
        boxVectors = self.topology.getPeriodicBoxVectors()
        if boxVectors is not None:
            sys.setDefaultPeriodicBoxVectors(*boxVectors)
        else:
            raise ValueError("periodicBoxVectors must be set")

        # custom non-bonded force to deal with martini's
        # use of epsilon_r
        es = mm.CustomNonbondedForce(
            f"step(rcut-r) * ES;"
            f"ES = f/epsilon_r*q1*q2 * (1/r + krf * r^2 - crf);"
            f"crf = 1 / rcut + krf * rcut^2;"
            f"krf = 1 / (2 * rcut^3);"
            f"epsilon_r = {epsilon_r};"
            f"f = 138.935458;"
            f"rcut={nonbondedCutoff.value_in_unit(unit.nanometers)};"
        )
        es.addPerParticleParameter("q")
        es.setNonbondedMethod(mm.CustomNonbondedForce.CutoffPeriodic)
        es.setCutoffDistance(nonbondedCutoff.value_in_unit(unit.nanometer))
        sys.addForce(es)

        # custom non-bonded force to add in the electrostatic
        # terms for self and excluded interactions
        elecSelfExclForce = mm.CustomBondForce(
            f"step(rcut-r) * ES;"
            f"ES = f*qprod/epsilon_r * (krf * r^2 - crf);"
            f"crf = 1 / rcut + krf * rcut^2;"
            f"krf = 1 / (2 * rcut^3);"
            f"epsilon_r = {epsilon_r};"
            f"f = 138.935458;"
            f"rcut={nonbondedCutoff.value_in_unit(unit.nanometers)};"
        )
        elecSelfExclForce.addPerBondParameter("qprod")
        sys.addForce(elecSelfExclForce)

        # custom LJ force using a lookup table to implment
        # nbfix-like terms
        lj = mm.CustomNonbondedForce(
            "step(rcut-r)*(energy - corr);"
            "energy = (C12(type1, type2) / r^12 - C6(type1, type2) / r^6);"
            "corr = (C12(type1, type2) / rcut^12 - C6(type1, type2) / rcut^6);"
            f"rcut={nonbondedCutoff.value_in_unit(unit.nanometers)};"
        )
        lj.addPerParticleParameter("type")
        lj.setNonbondedMethod(mm.CustomNonbondedForce.CutoffPeriodic)
        lj.setCutoffDistance(nonbondedCutoff.value_in_unit(unit.nanometer))
        sys.addForce(lj)

        # custom Lbonded force to handle exceptions
        ljExceptForce = mm.CustomBondForce(
            "step(rcut-r) * (energy - corr);"
            "energy = (C12/r^12 - C6/r^6);"
            "corr = (C12/rcut^12 - C6/rcut^6);"
            f"rcut={nonbondedCutoff.value_in_unit(unit.nanometers)};"
        )
        ljExceptForce.addPerBondParameter("C12")
        ljExceptForce.addPerBondParameter("C6")
        sys.addForce(ljExceptForce)

        allExceptions = []
        allCharges = []

        dihedralTypeTable, wildcardDihedralTypes = self._buildDihLookupTable()

        # build a lookup table mapping atom types into integer indices
        # that will later be used to lookup LJ combinations
        used_atom_types = set()
        for moleculeName, _ in self._molecules:
            moleculeType = self._moleculeTypes[moleculeName]
            for atom in moleculeType.atoms:
                used_atom_types.add(atom[1])
        atom_type_map = {k: i for i, k in enumerate(used_atom_types)}

        # now we need to setup our tables of C6 and C12
        n_types = len(atom_type_map)
        C6 = []
        C12 = []
        for i, type_i in enumerate(atom_type_map):
            for j, type_j in enumerate(atom_type_map):
                params = self._nonbondTypes[tuple(sorted([type_i, type_j]))]
                c6 = float(params[3])
                c12 = float(params[4])
                C6.append(c6)
                C12.append(c12)
        lj.addTabulatedFunction("C6", mm.Discrete2DFunction(n_types, n_types, C6))
        lj.addTabulatedFunction("C12", mm.Discrete2DFunction(n_types, n_types, C12))

        # Loop over molecules and create the specified number of each type.
        for moleculeName, moleculeCount in self._molecules:
            moleculeType = self._moleculeTypes[moleculeName]
            for i in range(moleculeCount):
                baseAtomIndex = sys.getNumParticles()
                atomTypes = [atom[1] for atom in moleculeType.atoms]
                try:
                    bondedTypes = [self._atomTypes[t][1] for t in atomTypes]
                except KeyError as e:
                    raise ValueError("Unknown atom type: " + e.message)
                bondedTypes = [
                    b if b is not None else a for a, b in zip(atomTypes, bondedTypes)
                ]

                # add bonded parameters
                self._addAtomsToSystem(sys, moleculeType)
                self._addBondsToSystem(sys, moleculeType, bondedTypes, baseAtomIndex)
                self._addAnglesToSystem(sys, moleculeType, bondedTypes, baseAtomIndex)
                self._addTorsionToSystem(
                    sys,
                    moleculeType,
                    bondedTypes,
                    dihedralTypeTable,
                    wildcardDihedralTypes,
                    baseAtomIndex,
                )
                self._addCmapToSystem(sys, moleculeType, bondedTypes, baseAtomIndex)

                exceptions, charges = self._setnonbondedParams(
                    es,
                    elecSelfExclForce,
                    lj,
                    moleculeType,
                    baseAtomIndex,
                    atomTypes,
                    atom_type_map,
                )
                allExceptions.extend(exceptions)
                allCharges.extend(charges)

                # Record virtual sites
                for fields in moleculeType.vsites2:
                    atoms = [int(x) - 1 for x in fields[:3]]
                    c1 = float(fields[4])
                    vsite = mm.TwoParticleAverageSite(
                        baseAtomIndex + atoms[1], baseAtomIndex + atoms[2], (1 - c1), c1
                    )
                    sys.setVirtualSite(baseAtomIndex + atoms[0], vsite)

                # Add explicitly specified constraints.
                for fields in moleculeType.constraints:
                    atoms = [int(x) - 1 for x in fields[:2]]
                    length = float(fields[2])
                    sys.addConstraint(
                        baseAtomIndex + atoms[0], baseAtomIndex + atoms[1], length
                    )

        if allExceptions:
            # build a map of unique exceptions
            # process in order, so that later entries trump earlier ones
            except_map = defaultdict(list)
            for exception in allExceptions:
                i, j, q, c6, c12 = exception
                if i < j:
                    except_map[(i, j)] = [q, c6, c12]
                else:
                    except_map[(j, i)] = [q, c6, C12]

            # add in all of the exclusions
            for i, j in except_map:
                es.addExclusion(i, j)
                lj.addExclusion(i, j)
                elecSelfExclForce.addBond(i, j, [allCharges[i] * allCharges[j]])
                ljExceptForce.addBond(i, j, [c6, c12])

        if removeCMMotion:
            sys.addForce(mm.CMMotionRemover())
        return sys


def _defaultGromacsIncludeDir():
    """Find the location where gromacs #include files are referenced from, by
    searching for (1) gromacs environment variables, (2) for the gromacs binary
    'pdb2gmx' or 'gmx' in the PATH, or (3) just using the default gromacs
    install location, /usr/local/gromacs/share/gromacs/top """
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
