"""Molecule data and atom-property preparation for DB2 generation."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TypeAlias


Coordinate: TypeAlias = tuple[float, float, float]
BondedAtom: TypeAlias = tuple[int, str]


@dataclass(slots=True)
class MoleculeData:
    """Topology and conformer coordinates consumed by the DB2 converter.

    Atom indices in ``atom_bonds`` are zero-based. The DB2-facing atom and
    bond numbers are stored separately because the output format is one-based.
    A topology can be prepared once and copied cheaply for several conformer
    batches with :meth:`copy_topology`.
    """

    name: str = "fake"
    protein_name: str = "fake"
    smiles: str = "fake"
    long_name: str = "fake"

    atom_numbers: list[int] = field(default_factory=list)
    atom_names: list[str] = field(default_factory=list)
    atom_types: list[str] = field(default_factory=list)
    atom_bonds: list[list[BondedAtom]] = field(default_factory=list)

    bond_numbers: list[int] = field(default_factory=list)
    bond_starts: list[int] = field(default_factory=list)
    bond_ends: list[int] = field(default_factory=list)
    bond_types: list[str] = field(default_factory=list)

    conformers: list[list[Coordinate]] = field(default_factory=list)
    input_energies: list[float] = field(default_factory=list)
    input_total_strain: list[float] = field(default_factory=list)
    input_max_strain: list[float] = field(default_factory=list)
    input_hydrogen_states: list[int] = field(default_factory=list)

    dock_atom_types: list[int] = field(default_factory=list)
    color_ids: list[int] = field(default_factory=list)

    def copy_topology(self) -> MoleculeData:
        """Return an independent topology without conformer-specific data."""

        return MoleculeData(
            name=self.name,
            protein_name=self.protein_name,
            smiles=self.smiles,
            long_name=self.long_name,
            atom_numbers=list(self.atom_numbers),
            atom_names=list(self.atom_names),
            atom_types=list(self.atom_types),
            atom_bonds=[list(neighbors) for neighbors in self.atom_bonds],
            bond_numbers=list(self.bond_numbers),
            bond_starts=list(self.bond_starts),
            bond_ends=list(self.bond_ends),
            bond_types=list(self.bond_types),
            dock_atom_types=list(self.dock_atom_types),
            color_ids=list(self.color_ids),
        )


DOCK_ATOM_TYPES: dict[str, int] = {
    "C.3": 5,
    "C.2": 1,
    "C.ar": 1,
    "C.1": 1,
    "N.3": 10,
    "N.2": 8,
    "N.1": 8,
    "O.3": 12,
    "O.2": 11,
    "S.3": 14,
    "N.ar": 8,
    "P.3": 13,
    "H": 6,
    "H-C": 7,
    "Br": 17,
    "Cl": 16,
    "F": 15,
    "I": 18,
    "S.2": 14,
    "N.pl3": 8,
    "LP": 25,
    "Na": 19,
    "K": 19,
    "Ca": 21,
    "Li": 20,
    "Al": 20,
    "Du": 25,
    "Du.C": 25,
    "Si": 24,
    "N.am": 8,
    "S.o": 14,
    "S.O": 14,
    "S.o2": 14,
    "S.O2": 14,
    "N.4": 9,
    "O.co2": 11,
    "C.cat": 1,
    "H.spc": 6,
    "O.spc": 11,
    "H.t3p": 6,
    "O.t3p": 11,
    "ANY": 25,
    "HEV": 25,
    "HET": 25,
    "HAL": 25,
    "Mg": 20,
    "Cr.oh": 25,
    "Cr.th": 25,
    "Se": 25,
    "Fe": 25,
    "Cu": 25,
    "Zn": 26,
    "Sn": 25,
    "Mo": 25,
    "Mn": 25,
    "Co.oh": 25,
}

COLOR_IDS: dict[str, int] = {
    "positive": 1,
    "negative": 2,
    "acceptor": 3,
    "donor": 4,
    "ester_o": 5,
    "amide_o": 6,
    "neutral": 7,
}

# Rules are evaluated in order and the final matching rule wins. A two-item
# rule matches an atom type; a four-item rule adds a graph-distance condition.
COLOR_RULES: tuple[tuple[object, ...], ...] = (
    ("N.4", "positive"),
    ("O.co2", "negative"),
    ("O.2", "acceptor"),
    ("O.3", "acceptor"),
    ("S.2", "acceptor"),
    ("N.ar", "acceptor"),
    ("O.3", 1, "P.3", "negative"),
    ("O.3", 1, "S.o2", "negative"),
    ("N.am", 1, "H", "donor"),
    ("N.pl3", 1, "H", "donor"),
    ("N.2", -1, "H", "acceptor"),
    ("N.2", -1, "C.3", "acceptor"),
    ("N.2", 1, "H", "donor"),
    ("N.ar", -1, "H", "acceptor"),
    ("N.ar", -1, "C.3", "acceptor"),
    ("N.ar", 1, "H", "donor"),
    ("O.3", 1, "H", "donor"),
    ("O.2", 2, "O.3", "ester_o"),
    ("O.2", 2, "N.pl3", "amide_o"),
    ("O.2", 2, "N.am", "amide_o"),
    ("O.2", 2, "N.3", "amide_o"),
)


def prepare_molecule_for_db2(molecule: MoleculeData) -> MoleculeData:
    """Complete per-conformer metadata and calculate DB2 atom properties.

    The metadata lists and the calculated ``dock_atom_types``/``color_ids``
    fields are updated in place, matching the historical entry-point side
    effects. Graph neighborhoods are calculated once and reused by every color
    rule instead of repeatedly traversing the molecular graph.
    """

    conformer_count = len(molecule.conformers)
    while len(molecule.input_energies) < conformer_count:
        molecule.input_energies.append(9999.99)
        molecule.input_total_strain.append(0.0)
        molecule.input_max_strain.append(0.0)
        molecule.input_hydrogen_states.append(0)

    graph_levels = _build_graph_levels(molecule.atom_bonds, maximum_depth=2)
    molecule.dock_atom_types = _assign_dock_atom_types(molecule, graph_levels)
    molecule.color_ids = _assign_color_ids(molecule, graph_levels)
    return molecule


def _build_graph_levels(
    atom_bonds: list[list[BondedAtom]], maximum_depth: int
) -> list[list[list[int]]]:
    """Return atoms at each exact shortest-path depth for every atom."""

    levels_by_atom: list[list[list[int]]] = []
    for start_atom in range(len(atom_bonds)):
        levels = [[start_atom]]
        visited = {start_atom}
        frontier = [start_atom]
        for _ in range(maximum_depth):
            next_frontier: list[int] = []
            for atom_index in frontier:
                for neighbor_index, _bond_type in atom_bonds[atom_index]:
                    if neighbor_index not in visited:
                        visited.add(neighbor_index)
                        next_frontier.append(neighbor_index)
            levels.append(next_frontier)
            frontier = next_frontier
        levels_by_atom.append(levels)
    return levels_by_atom


def _has_bonded_type(
    atom_types: list[str],
    graph_levels: list[list[list[int]]],
    atom_index: int,
    type_fragment: str,
    bonds_away: int,
) -> bool:
    return any(
        type_fragment in atom_types[other_index]
        for other_index in graph_levels[atom_index][bonds_away]
    )


def _assign_dock_atom_types(
    molecule: MoleculeData, graph_levels: list[list[list[int]]]
) -> list[int]:
    dock_types: list[int] = []
    for atom_index, atom_type in enumerate(molecule.atom_types):
        if atom_type == "H" and _has_bonded_type(
            molecule.atom_types, graph_levels, atom_index, "C", 1
        ):
            dock_types.append(DOCK_ATOM_TYPES["H-C"])
        else:
            dock_types.append(DOCK_ATOM_TYPES[atom_type])
    return dock_types


def _assign_color_ids(
    molecule: MoleculeData, graph_levels: list[list[list[int]]]
) -> list[int]:
    color_ids: list[int] = []
    for atom_index, atom_type in enumerate(molecule.atom_types):
        color_name = "neutral"
        for rule in COLOR_RULES:
            if not atom_type.startswith(rule[0]):
                continue
            if len(rule) == 2:
                color_name = rule[1]
                continue

            bonds_away = rule[1]
            bonded_type = rule[2]
            matches_bond = _has_bonded_type(
                molecule.atom_types,
                graph_levels,
                atom_index,
                bonded_type,
                abs(bonds_away),
            )
            if (bonds_away == -1 and not matches_bond) or (
                bonds_away > 0 and matches_bond
            ):
                color_name = rule[3]

        if atom_type == "O.3":
            bonded_to_phosphorus_or_sulfur = _has_bonded_type(
                molecule.atom_types, graph_levels, atom_index, "P.3", 1
            ) or _has_bonded_type(
                molecule.atom_types, graph_levels, atom_index, "S.o2", 1
            )
            bonded_to_carbon = _has_bonded_type(
                molecule.atom_types, graph_levels, atom_index, "C.3", 1
            ) or _has_bonded_type(
                molecule.atom_types, graph_levels, atom_index, "C.2", 1
            )
            if bonded_to_phosphorus_or_sulfur and bonded_to_carbon:
                color_name = "acceptor"

        color_ids.append(COLOR_IDS[color_name])
    return color_ids
