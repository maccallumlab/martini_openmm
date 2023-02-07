from collections import defaultdict
from typing import List, Dict, Union
import simtk.unit as unit


class LinearSite:
    atom_weights: Dict[int, float]

    def __init__(self, atom_weight_dict):
        self.atom_weights = atom_weight_dict

    def __repr__(self):
        return f"LinearSite({self.atom_weights})"


class COMLinearSite:
    atoms: List[int]

    def __init__(self, atoms):
        self.atoms = atoms

    def to_linear(self, system, offset):
        masses = [
            system.getParticleMass(i + offset).value_in_unit(unit.dalton)
            for i in self.atoms
        ]
        total = sum(masses)
        weights = [m / total for m in masses]
        site_dict = {atom: weight for atom, weight in zip(self.atoms, weights)}
        site = LinearSite(site_dict)
        return site

    def __repr__(self):
        return f"COMLinearSite({self.atoms})"


class NonLinearSite:
    pass


class OutOfPlane(NonLinearSite):
    atom1: int
    atom2: int
    atom3: int
    a: float
    b: float
    c: float

    def __init__(self, atom1, atom2, atom3, a, b, c):
        self.atom1 = atom1
        self.atom2 = atom2
        self.atom3 = atom3
        self.a = a
        self.b = b
        self.c = c

    def __repr__(self):
        return f"OutOfPlane({self.atom1}, {self.atom2}, {self.atom3}, {self.a}, {self.b}, {self.c})"


class NormalizedInPlaneSite(NonLinearSite):
    atom1: int
    atom2: int
    atom3: int
    a: float
    d: float

    def __init__(self, atom1, atom2, atom3, a, d) -> None:
        self.atom1 = atom1
        self.atom2 = atom2
        self.atom3 = atom3
        self.a = a
        self.d = d

    def __repr__(self):
        return f"NormalizedInPlaneSite({self.atom1}, {self.atom2}, {self.atom3}, {self.a}, {self.d})"


class NormalizedInPlaneTwoParticleSite(NonLinearSite):
    atom1: int
    atom2: int
    a: float

    def __init__(self, atom1, atom2, a) -> None:
        self.atom1 = atom1
        self.atom2 = atom2
        self.a = a

    def __repr__(self):
        return f"NormalizedInPlaneTwoParticleSite({self.atom1}, {self.atom2}, {self.a})"


class VSiteManager:
    vsites: Dict[int, Union[LinearSite, NonLinearSite]]

    def __init__(self):
        self.vsites = {}

    def add(self, index, site):
        if index in self.vsites:
            raise ValueError(f"Tried to add more than one vsite for particle {index}.")
        self.vsites[index] = site

    def convert_com_to_linear(self, system, offset):
        for index, site in self.vsites.items():
            if isinstance(site, COMLinearSite):
                self.vsites[index] = site.to_linear(system, offset)

    def iter(self):
        for index, site in self.vsites.items():
            if isinstance(site, NonLinearSite):
                yield index, site
            else:
                site = self.flatten_site(site)
                yield index, site

    def flatten_site(self, site):
        if isinstance(site, NonLinearSite):
            raise RuntimeError("A linear vsite cannot depend on a non-linear vsite.")

        # If this is a COM site, find the masses and convert
        # to a LinearSite
        if isinstance(site, COMLinearSite):
            raise RuntimeError(
                "All COM sites should have been converted to LinearSites."
            )

        from_atoms = defaultdict(list)

        for atom, weight in site.atom_weights.items():
            # Atom is itself a vsite.
            # We recrursively flatten and multiply by the weight.
            if atom in self.vsites:
                flattened = self.flatten_site(self.vsites[atom])
                for f, w in flattened.atom_weights.items():
                    from_atoms[f].append(weight * w)

            # atom is just a regular atom
            else:
                from_atoms[atom].append(weight)

        # Now need to sum the coefficients for each from_atom.
        from_atoms = {atom: sum(weights) for atom, weights in from_atoms.items()}
        return LinearSite(from_atoms)
