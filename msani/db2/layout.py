"""Build and serialize the coordinate layout used by the DB2 format."""

from __future__ import annotations

from dataclasses import dataclass
from itertools import product
import math
from typing import Any

from msani.db2.molecule import Coordinate, MoleculeData


SET_CONFORMERS_PER_LINE = 8

_MOORE_NEIGHBORHOOD = tuple(product(range(-1, 2), repeat=3))
_MOORE_NEIGHBORHOOD = tuple(
    offset for offset in _MOORE_NEIGHBORHOOD if offset != (0, 0, 0)
)
_EXTENDED_MOORE_FACES = (
    tuple((-2, second, third) for second, third in product(range(-1, 2), repeat=2)),
    tuple((2, second, third) for second, third in product(range(-1, 2), repeat=2)),
    tuple((first, -2, third) for first, third in product(range(-1, 2), repeat=2)),
    tuple((first, 2, third) for first, third in product(range(-1, 2), repeat=2)),
    tuple((first, second, -2) for first, second in product(range(-1, 2), repeat=2)),
    tuple((first, second, 2) for first, second in product(range(-1, 2), repeat=2)),
)


@dataclass(frozen=True, slots=True)
class Db2Coordinate:
    """One unique coordinate emitted by the DB2 ``X`` section."""

    atom_index: int
    conformation_number: int
    xyz: Coordinate


@dataclass(frozen=True, slots=True)
class Db2Layout:
    """DB2-specific coordinate, conformation, and input-set relationships."""

    position_counts: tuple[int, ...]
    rigid_structure_ids: tuple[int, ...]
    rigid_atom_indices: tuple[int, ...]
    heavy_rigid_atom_indices: tuple[int, ...]
    coordinates: tuple[Db2Coordinate, ...]
    conformer_coordinate_ranges: tuple[tuple[int, int], ...]
    conformers_by_input_set: dict[int, tuple[int, ...]]

    @property
    def conformation_count(self) -> int:
        return len(self.conformer_coordinate_ranges)


@dataclass(slots=True)
class _PositionGroup:
    """Atoms sharing the same membership across input conformer sets."""

    atom_indices: list[int]
    representative_coordinates: list[Coordinate]


@dataclass(slots=True)
class _SpatialBucket:
    conformer_indices: list[int]
    visited: bool = False


class _DisjointSets:
    """Small union/find implementation for fixed molecular components."""

    def __init__(self) -> None:
        self._parents: dict[int, int] = {}
        self._ranks: dict[int, int] = {}

    def find(self, item: int) -> int:
        if item not in self._parents:
            self._parents[item] = item
            self._ranks[item] = 0
            return item

        path = [item]
        parent = self._parents[item]
        while parent != path[-1]:
            path.append(parent)
            parent = self._parents[parent]
        for path_item in path[:-1]:
            self._parents[path_item] = parent
        return parent

    def union(self, first: int, second: int) -> None:
        first_parent = self.find(first)
        second_parent = self.find(second)
        if first_parent == second_parent:
            return
        if self._ranks[first_parent] < self._ranks[second_parent]:
            self._parents[first_parent] = second_parent
            return
        self._parents[second_parent] = first_parent
        if self._ranks[first_parent] == self._ranks[second_parent]:
            self._ranks[first_parent] += 1

    def groups(self) -> list[list[int]]:
        groups_by_parent: dict[int, list[int]] = {}
        for item in self._parents:
            self.find(item)
        for item, parent in self._parents.items():
            groups_by_parent.setdefault(parent, []).append(item)
        return list(groups_by_parent.values())


class _SpatialPositionClusterer:
    """Cluster conformer positions for one atom within a distance tolerance."""

    def __init__(self, tolerance: float) -> None:
        self.tolerance_squared = tolerance**2
        self.bucket_size = tolerance / math.sqrt(3)
        self.extra_width = tolerance - self.bucket_size

    def cluster(
        self, positions: list[Coordinate]
    ) -> list[tuple[tuple[int, ...], Coordinate]]:
        buckets: dict[tuple[int, int, int], _SpatialBucket] = {}
        for conformer_index, xyz in enumerate(positions):
            bucket_key = (
                math.floor(xyz[0] / self.bucket_size),
                math.floor(xyz[1] / self.bucket_size),
                math.floor(xyz[2] / self.bucket_size),
            )
            bucket = buckets.get(bucket_key)
            if bucket is None:
                buckets[bucket_key] = _SpatialBucket([conformer_index])
            else:
                bucket.conformer_indices.append(conformer_index)

        clusters: list[tuple[tuple[int, ...], Coordinate]] = []
        for bucket_key, bucket in buckets.items():
            if not bucket.conformer_indices:
                continue

            representative = positions[bucket.conformer_indices[0]]
            cluster_members = list(bucket.conformer_indices)

            def absorb_nearby(neighbor_key: tuple[int, int, int]) -> None:
                neighbor_bucket = buckets.get(neighbor_key)
                if neighbor_bucket is None or neighbor_bucket.visited:
                    return

                remaining_indices: list[int] = []
                for neighbor_index in neighbor_bucket.conformer_indices:
                    neighbor = positions[neighbor_index]
                    x_difference = representative[0] - neighbor[0]
                    y_difference = representative[1] - neighbor[1]
                    z_difference = representative[2] - neighbor[2]
                    if (
                        x_difference * x_difference
                        + y_difference * y_difference
                        + z_difference * z_difference
                        <= self.tolerance_squared
                    ):
                        cluster_members.append(neighbor_index)
                    else:
                        remaining_indices.append(neighbor_index)
                neighbor_bucket.conformer_indices = remaining_indices

            for offset in _MOORE_NEIGHBORHOOD:
                absorb_nearby(
                    (
                        bucket_key[0] + offset[0],
                        bucket_key[1] + offset[1],
                        bucket_key[2] + offset[2],
                    )
                )
            for face_index in self._faces_near_position(representative, bucket_key):
                for offset in _EXTENDED_MOORE_FACES[face_index]:
                    absorb_nearby(
                        (
                            bucket_key[0] + offset[0],
                            bucket_key[1] + offset[1],
                            bucket_key[2] + offset[2],
                        )
                    )

            bucket.visited = True
            clusters.append((tuple(sorted(cluster_members)), representative))
        return clusters

    def _faces_near_position(
        self, xyz: Coordinate, bucket_key: tuple[int, int, int]
    ) -> list[int]:
        x_offset = xyz[0] - bucket_key[0] * self.bucket_size
        y_offset = xyz[1] - bucket_key[1] * self.bucket_size
        z_offset = xyz[2] - bucket_key[2] * self.bucket_size
        faces: list[int] = []
        if x_offset < self.extra_width:
            faces.append(0)
        if x_offset > self.bucket_size - self.extra_width:
            faces.append(1)
        if y_offset < self.extra_width:
            faces.append(2)
        if y_offset > self.bucket_size - self.extra_width:
            faces.append(3)
        if z_offset < self.extra_width:
            faces.append(4)
        if z_offset > self.bucket_size - self.extra_width:
            faces.append(5)
        return faces


def build_db2_layout(
    molecule: MoleculeData, tolerance: float = 0.001
) -> Db2Layout:
    """Build the DB2 coordinate layout without changing coordinates.

    Positions within ``tolerance`` are represented by the first coordinate in
    their historical spatial cluster. Rigid structures describe atoms that are
    joined by non-rotatable bonds or belong to intersecting graph cycles.
    """

    if not molecule.conformers:
        raise ValueError("DB2 generation requires at least one conformer")

    rigid_structure_ids = _assign_rigid_structure_ids(
        len(molecule.atom_numbers), molecule.atom_bonds
    )
    position_counts, position_groups = _cluster_atom_positions(
        molecule.conformers, tolerance
    )
    rigid_atom_indices = _find_largest_fixed_component(
        position_counts, molecule.atom_bonds
    )
    coordinates, coordinate_ranges, conformers_by_set = _assemble_layout_records(
        position_groups,
        rigid_structure_ids,
        input_set_count=len(molecule.conformers),
    )

    rigid_atom_set = set(rigid_atom_indices)
    heavy_rigid_atom_indices = tuple(
        atom_index
        for atom_index in rigid_atom_set
        if "H" not in molecule.atom_types[atom_index]
    )
    return Db2Layout(
        position_counts=tuple(position_counts),
        rigid_structure_ids=tuple(rigid_structure_ids),
        rigid_atom_indices=tuple(rigid_atom_indices),
        heavy_rigid_atom_indices=heavy_rigid_atom_indices,
        coordinates=tuple(coordinates),
        conformer_coordinate_ranges=tuple(coordinate_ranges),
        conformers_by_input_set={
            set_index: tuple(conformer_numbers)
            for set_index, conformer_numbers in conformers_by_set.items()
        },
    )


def _assign_rigid_structure_ids(
    atom_count: int, atom_bonds: list[list[tuple[int, str]]]
) -> list[int]:
    """Map atoms to rigid segments using cycles and non-single bonds."""

    if atom_count == 0:
        raise ValueError("DB2 generation requires at least one atom")

    cycles_by_atom: list[list[int]] = [[] for _ in range(atom_count)]
    parent = [0] * atom_count
    visit_state = [0] * atom_count
    cycle_count = 0

    def find_graph_cycles(atom_index: int, previous_atom: int) -> None:
        nonlocal cycle_count
        if visit_state[atom_index] == 2:
            return
        if visit_state[atom_index] == 1:
            current_atom = previous_atom
            while current_atom != atom_index:
                cycles_by_atom[current_atom].append(cycle_count)
                current_atom = parent[current_atom]
            cycles_by_atom[current_atom].append(cycle_count)
            cycle_count += 1
            return

        visit_state[atom_index] = 1
        parent[atom_index] = previous_atom
        for neighbor_index, _bond_type in atom_bonds[atom_index]:
            if neighbor_index != previous_atom:
                find_graph_cycles(neighbor_index, atom_index)
        visit_state[atom_index] = 2

    # The historical implementation started only at atom zero. Visiting every
    # component retains connected-molecule behavior and makes disconnected
    # topologies deterministic instead of leaving them partially initialized.
    for atom_index in range(atom_count):
        if visit_state[atom_index] == 0:
            find_graph_cycles(atom_index, -1)

    for atom_index, neighbors in enumerate(atom_bonds):
        if len(neighbors) == 1:
            neighbor_index = neighbors[0][0]
            cycles_by_atom[atom_index].append(cycle_count)
            cycles_by_atom[neighbor_index].append(cycle_count)
            cycle_count += 1
            continue

        for neighbor_index, bond_type in neighbors:
            if neighbor_index < atom_index:
                continue
            shares_cycle = bool(
                set(cycles_by_atom[atom_index]).intersection(
                    cycles_by_atom[neighbor_index]
                )
            )
            if not shares_cycle and bond_type != "1":
                cycles_by_atom[atom_index].append(cycle_count)
                cycles_by_atom[neighbor_index].append(cycle_count)
                cycle_count += 1

    intersecting_cycles: list[set[int]] = [set() for _ in range(cycle_count)]
    for atom_cycles in cycles_by_atom:
        for first_offset, first_cycle in enumerate(atom_cycles):
            for second_cycle in atom_cycles[first_offset + 1 :]:
                intersecting_cycles[first_cycle].add(second_cycle)
                intersecting_cycles[second_cycle].add(first_cycle)

    rigid_id_by_cycle = [0] * cycle_count
    visited_cycles = [False] * cycle_count
    rigid_structure_count = 0

    def merge_intersecting_cycles(cycle_index: int) -> None:
        if visited_cycles[cycle_index]:
            return
        visited_cycles[cycle_index] = True
        rigid_id_by_cycle[cycle_index] = rigid_structure_count
        for intersecting_cycle in intersecting_cycles[cycle_index]:
            merge_intersecting_cycles(intersecting_cycle)

    for cycle_index in range(cycle_count):
        if not visited_cycles[cycle_index]:
            merge_intersecting_cycles(cycle_index)
            rigid_structure_count += 1

    rigid_structure_ids = [0] * atom_count
    for atom_index, atom_cycles in enumerate(cycles_by_atom):
        if atom_cycles:
            rigid_structure_ids[atom_index] = rigid_id_by_cycle[atom_cycles[0]]
        else:
            rigid_structure_ids[atom_index] = rigid_structure_count
            rigid_structure_count += 1
    return rigid_structure_ids


def _cluster_atom_positions(
    conformers: list[list[Coordinate]], tolerance: float
) -> tuple[list[int], dict[tuple[int, ...], _PositionGroup]]:
    input_set_count = len(conformers)
    atom_count = len(conformers[0])
    all_input_sets = tuple(range(input_set_count))
    tolerance_squared = tolerance * tolerance
    clusterer = _SpatialPositionClusterer(tolerance)
    position_counts: list[int] = []
    position_groups: dict[tuple[int, ...], _PositionGroup] = {}

    for atom_index in range(atom_count):
        atom_positions = [conformer[atom_index] for conformer in conformers]
        representative = atom_positions[0]
        fixed = True
        for xyz in atom_positions[1:]:
            x_difference = representative[0] - xyz[0]
            y_difference = representative[1] - xyz[1]
            z_difference = representative[2] - xyz[2]
            if (
                x_difference * x_difference
                + y_difference * y_difference
                + z_difference * z_difference
                > tolerance_squared
            ):
                fixed = False
                break

        if fixed:
            group = position_groups.get(all_input_sets)
            if group is None:
                position_groups[all_input_sets] = _PositionGroup(
                    [atom_index], [representative]
                )
            else:
                group.atom_indices.append(atom_index)
                group.representative_coordinates.append(representative)
            position_counts.append(1)
            continue

        spatial_clusters = clusterer.cluster(atom_positions)
        for input_sets, cluster_representative in spatial_clusters:
            group = position_groups.get(input_sets)
            if group is None:
                position_groups[input_sets] = _PositionGroup(
                    [atom_index], [cluster_representative]
                )
            else:
                group.atom_indices.append(atom_index)
                group.representative_coordinates.append(cluster_representative)
        position_counts.append(len(spatial_clusters))

    return position_counts, position_groups


def _find_largest_fixed_component(
    position_counts: list[int], atom_bonds: list[list[tuple[int, str]]]
) -> list[int]:
    components = _DisjointSets()
    for atom_index, position_count in enumerate(position_counts):
        if position_count != 1:
            continue
        for neighbor_index, _bond_type in atom_bonds[atom_index]:
            if position_counts[neighbor_index] == 1:
                components.union(atom_index, neighbor_index)

    largest_component: list[int] | None = None
    for component in components.groups():
        if largest_component is None or len(component) > len(largest_component):
            largest_component = component
    if largest_component is None:
        raise ValueError(
            "DB2 layout requires at least two bonded atoms fixed across conformers"
        )
    return largest_component


def _assemble_layout_records(
    position_groups: dict[tuple[int, ...], _PositionGroup],
    rigid_structure_ids: list[int],
    input_set_count: int,
) -> tuple[
    list[Db2Coordinate],
    list[tuple[int, int]],
    dict[int, list[int]],
]:
    coordinates: list[Db2Coordinate] = []
    coordinate_ranges: list[tuple[int, int]] = []
    conformers_by_set = {set_index: [] for set_index in range(input_set_count)}
    conformation_number = 0

    groups_by_decreasing_membership = sorted(
        position_groups.items(), key=lambda item: -len(item[0])
    )
    for input_sets, group in groups_by_decreasing_membership:
        # Sorting rigid IDs determines conformation boundaries. Atom and
        # coordinate record order intentionally remains the historical atom
        # processing order so DB2 output stays byte-for-byte compatible.
        sorted_rigid_ids = sorted(
            rigid_structure_ids[atom_index] for atom_index in group.atom_indices
        )
        previous_rigid_id: int | None = None
        range_start = 0
        for group_offset, rigid_id in enumerate(sorted_rigid_ids):
            if rigid_id != previous_rigid_id:
                if previous_rigid_id is not None:
                    coordinate_ranges.append((range_start, len(coordinates) - 1))
                range_start = len(coordinates)
                conformation_number += 1
                for set_index in input_sets:
                    conformers_by_set[set_index].append(conformation_number - 1)

            coordinates.append(
                Db2Coordinate(
                    atom_index=group.atom_indices[group_offset],
                    conformation_number=conformation_number,
                    xyz=group.representative_coordinates[group_offset],
                )
            )
            previous_rigid_id = rigid_id

        coordinate_ranges.append((range_start, len(coordinates) - 1))

    return coordinates, coordinate_ranges, conformers_by_set


def serialize_db2(
    molecule: MoleculeData,
    solvation_data: Any,
    layout: Db2Layout,
) -> str:
    """Serialize a prepared molecule and its DB2 layout."""

    lines: list[str] = []
    lines.append(
        "M %16s %9s %3d %3d %6d %6d %6d %6d %6d %6d\n"
        % (
            molecule.name[-16:],
            molecule.protein_name[-9:],
            len(molecule.atom_numbers),
            len(molecule.bond_starts),
            len(layout.coordinates),
            layout.conformation_count,
            len(layout.conformers_by_input_set),
            len(layout.heavy_rigid_atom_indices),
            5,
            0,
        )
    )
    lines.append(
        "M %+9.4f %+10.3f %+10.3f %+10.3f %9.3f\n"
        % (
            solvation_data.totalCharge,
            solvation_data.totalPolarSolv,
            solvation_data.totalApolarSolv,
            solvation_data.totalSolv,
            solvation_data.totalSurface,
        )
    )
    lines.append("M %-76s\n" % molecule.smiles[-76:])
    lines.append("M %-76s\n" % molecule.long_name[-76:])
    lines.append("M %+10.4f\n" % 999.999)

    for atom_index, atom_number in enumerate(molecule.atom_numbers):
        lines.append(
            "A %3d %-4s %-5s %2d %2d %+9.4f %+10.3f %+10.3f %+10.3f %9.3f\n"
            % (
                atom_number,
                molecule.atom_names[atom_index],
                molecule.atom_types[atom_index],
                molecule.dock_atom_types[atom_index],
                molecule.color_ids[atom_index],
                solvation_data.charge[atom_index],
                solvation_data.polarSolv[atom_index],
                solvation_data.apolarSolv[atom_index],
                solvation_data.solv[atom_index],
                solvation_data.surface[atom_index],
            )
        )
    for bond_index, bond_number in enumerate(molecule.bond_numbers):
        lines.append(
            "B %3d %3d %3d %-2s\n"
            % (
                bond_number,
                molecule.bond_starts[bond_index],
                molecule.bond_ends[bond_index],
                molecule.bond_types[bond_index],
            )
        )
    for coordinate_number, coordinate in enumerate(layout.coordinates, 1):
        lines.append(
            "X %9d %3d %6d %+9.4f %+9.4f %+9.4f\n"
            % (
                coordinate_number,
                coordinate.atom_index + 1,
                coordinate.conformation_number,
                coordinate.xyz[0],
                coordinate.xyz[1],
                coordinate.xyz[2],
            )
        )
    for rigid_number, atom_index in enumerate(
        layout.heavy_rigid_atom_indices, 1
    ):
        xyz = molecule.conformers[0][atom_index]
        lines.append(
            "R %6d %2d %+9.4f %+9.4f %+9.4f\n"
            % (
                rigid_number,
                molecule.color_ids[atom_index],
                xyz[0],
                xyz[1],
                xyz[2],
            )
        )
    for conformation_number, (start, end) in enumerate(
        layout.conformer_coordinate_ranges, 1
    ):
        lines.append(
            "C %6d %9d %9d\n"
            % (conformation_number, start + 1, end + 1)
        )

    for output_set_number, input_set_index in enumerate(
        sorted(layout.conformers_by_input_set), 1
    ):
        conformer_numbers = [
            number + 1
            for number in layout.conformers_by_input_set[input_set_index]
        ]
        total_conformers = len(conformer_numbers)
        total_lines = math.ceil(total_conformers / SET_CONFORMERS_PER_LINE)
        lines.append(
            "S %6d %6d %3d %1d %1d %+11.3f %+11.3f\n"
            % (
                output_set_number,
                total_lines,
                total_conformers,
                0,
                _value_or_default(molecule.input_hydrogen_states, input_set_index, 0),
                _value_or_default(molecule.input_total_strain, input_set_index, 0.0),
                _value_or_default(molecule.input_max_strain, input_set_index, 0.0),
            )
        )
        for line_number in range(total_lines):
            line_conformers = conformer_numbers[
                line_number
                * SET_CONFORMERS_PER_LINE : (line_number + 1)
                * SET_CONFORMERS_PER_LINE
            ]
            line = "S %6d %6d %1d" % (
                output_set_number,
                line_number + 1,
                len(line_conformers),
            )
            line += "".join(" %6d" % number for number in line_conformers)
            lines.append(line + "\n")

    lines.append("E\n")
    return "".join(lines)


def _value_or_default(values: list[Any], index: int, default: Any) -> Any:
    """Return sparse optional conformer metadata without materializing defaults."""

    return values[index] if index < len(values) else default
