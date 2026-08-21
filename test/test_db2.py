"""Characterization and stage-level tests for MOL2-to-DB2 conversion."""

from __future__ import annotations

from hashlib import sha256
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest

from msani.db2 import mol2db2 as mol2db2_module
from msani.db2.hierarchy import build_conformer_hierarchy, serialize_db2
from msani.db2.molecule import MoleculeData, prepare_molecule_for_db2


def _molecule(
    name,
    atom_types,
    bonds,
    conformers,
    total_strain=None,
    max_strain=None,
    hydrogen_states=None,
):
    atom_count = len(atom_types)
    conformer_count = len(conformers)
    molecule = MoleculeData(
        name=name,
        protein_name="none",
        smiles=f"{name}_smiles",
        long_name=f"{name} long name",
        atom_numbers=list(range(1, atom_count + 1)),
        atom_names=[f"A{number}" for number in range(1, atom_count + 1)],
        atom_types=list(atom_types),
        atom_bonds=[[] for _ in range(atom_count)],
        conformers=[list(conformer) for conformer in conformers],
        input_energies=[9999.99] * conformer_count,
        input_total_strain=list(total_strain or [0.0] * conformer_count),
        input_max_strain=list(max_strain or [0.0] * conformer_count),
        input_hydrogen_states=list(hydrogen_states or [0] * conformer_count),
    )
    for bond_number, (start, end, bond_type) in enumerate(bonds, 1):
        molecule.bond_numbers.append(bond_number)
        molecule.bond_starts.append(start)
        molecule.bond_ends.append(end)
        molecule.bond_types.append(bond_type)
        molecule.atom_bonds[start - 1].append((end - 1, bond_type))
        molecule.atom_bonds[end - 1].append((start - 1, bond_type))
    return molecule


def _solvation(atom_count, charges=None):
    charges = list(charges or [0.0] * atom_count)
    polar = [-(index + 1) * 0.125 for index in range(atom_count)]
    apolar = [(index + 1) * 0.075 for index in range(atom_count)]
    surface = [(index + 1) * 1.25 for index in range(atom_count)]
    return SimpleNamespace(
        charge=charges,
        polarSolv=polar,
        apolarSolv=apolar,
        solv=[polar_value + apolar_value for polar_value, apolar_value in zip(polar, apolar)],
        surface=surface,
        totalCharge=sum(charges),
        totalPolarSolv=sum(polar),
        totalApolarSolv=sum(apolar),
        totalSolv=sum(polar) + sum(apolar),
        totalSurface=sum(surface),
    )


def _flexible_molecule():
    return _molecule(
        "flexible",
        ["C.3", "C.3", "O.3", "H"],
        [(1, 2, "1"), (2, 3, "1"), (3, 4, "1")],
        [
            [(0.0, 0.0, 0.0), (1.5, 0.0, 0.0), (2.5, 1.0, 0.0), (3.2, 1.4, 0.0)],
            [(0.0, 0.0, 0.0), (1.5, 0.0, 0.0), (2.5, -1.0, 0.0), (3.2, -1.4, 0.0)],
            [(0.0, 0.0, 0.0), (1.5, 0.0, 0.0), (2.5, 1.0, 0.0), (3.2, 1.4, 0.0)],
        ],
        total_strain=[0.0, 2.5, 0.25],
        max_strain=[0.0, 1.5, 0.1],
        hydrogen_states=[0, 2, 1],
    )


class TestDb2Conversion(unittest.TestCase):
    def test_supported_module_import_and_rigid_charged_output(self):
        molecule = _molecule(
            "rigid_charged",
            ["C.2", "O.co2", "O.co2"],
            [(1, 2, "2"), (1, 3, "1")],
            [[(0.0, 0.0, 0.0), (1.2, 0.0, 0.0), (-0.6, 1.0392304845, 0.0)]],
            total_strain=[1.25],
            max_strain=[0.75],
        )

        output = mol2db2_module.mol2db2(
            molecule, _solvation(3, [0.5, -0.75, -0.75])
        )

        self.assertEqual(
            sha256(output.encode()).hexdigest(),
            "452079221f3cb2a66d36eaf7403ab77ebde950ebf897dcab16e12e0e283c6c4f",
        )
        self.assertEqual(molecule.dock_atom_types, [1, 11, 11])
        self.assertEqual(molecule.color_ids, [7, 2, 2])
        self.assertTrue(output.endswith("E\n"))

    def test_flexible_hierarchy_stages_and_output(self):
        prepared = prepare_molecule_for_db2(_flexible_molecule())
        hierarchy = build_conformer_hierarchy(prepared)

        self.assertEqual(hierarchy.position_counts, (1, 1, 2, 2))
        self.assertEqual(hierarchy.rigid_atom_indices, (0, 1))
        self.assertEqual(hierarchy.heavy_rigid_atom_indices, (0, 1))
        self.assertEqual(
            hierarchy.conformer_coordinate_ranges,
            ((0, 1), (2, 3), (4, 5)),
        )
        self.assertEqual(
            hierarchy.conformers_by_input_set,
            {0: (0, 1), 1: (0, 2), 2: (0, 1)},
        )

        output = serialize_db2(prepared, _solvation(4), hierarchy)
        self.assertEqual(
            sha256(output.encode()).hexdigest(),
            "d3eb3256e8a615484ede97a088cabbe015e79a8e2e3adb493cb65fae743ad7d2",
        )
        self.assertEqual(sum(line.startswith("X ") for line in output.splitlines()), 6)
        self.assertEqual(sum(line.startswith("C ") for line in output.splitlines()), 3)

    def test_tolerance_clustering_matches_characterized_boundary_case(self):
        molecule = _molecule(
            "tolerance",
            ["C.3"] * 4,
            [(1, 2, "1"), (2, 3, "1"), (3, 4, "1")],
            [
                [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (2.0, 0.0, 0.0), (3.0, 0.0, 0.0)],
                [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (2.0005, 0.0, 0.0), (3.001, 0.0, 0.0)],
                [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (2.002, 0.0, 0.0), (3.003, 0.0, 0.0)],
            ],
        )

        output = mol2db2_module.mol2db2(molecule, _solvation(4), disttol=0.001)

        self.assertEqual(
            sha256(output.encode()).hexdigest(),
            "b7fea3065440afff58d489b67f79dbd0c60b1c36cf2f51645183a6a603a513df",
        )

    def test_missing_conformer_metadata_is_padded_in_place(self):
        molecule = _flexible_molecule()
        molecule.input_energies.clear()
        molecule.input_total_strain.clear()
        molecule.input_max_strain.clear()
        molecule.input_hydrogen_states.clear()

        mol2db2_module.mol2db2(molecule, _solvation(4))

        self.assertEqual(molecule.input_energies, [9999.99] * 3)
        self.assertEqual(molecule.input_total_strain, [0.0] * 3)
        self.assertEqual(molecule.input_max_strain, [0.0] * 3)
        self.assertEqual(molecule.input_hydrogen_states, [0] * 3)

    def test_supplied_clash_file_is_still_parsed_but_does_not_filter(self):
        molecule = _flexible_molecule()
        solution = _solvation(4)
        expected = mol2db2_module.mol2db2(molecule, solution)
        with TemporaryDirectory() as temporary_directory:
            clash_file = Path(temporary_directory) / "clash.txt"
            clash_file.write_text("min 2 1 H H 1.6\n")
            observed = mol2db2_module.mol2db2(
                molecule, solution, clashfile=clash_file
            )
        self.assertEqual(observed, expected)

    def test_unanchored_conformers_report_the_existing_unsupported_case(self):
        molecule = _molecule(
            "unanchored",
            ["C.3", "C.3"],
            [(1, 2, "1")],
            [
                [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0)],
                [(1.0, 0.0, 0.0), (2.0, 0.0, 0.0)],
            ],
        )
        with self.assertRaisesRegex(ValueError, "two bonded atoms fixed"):
            mol2db2_module.mol2db2(molecule, _solvation(2))


if __name__ == "__main__":
    unittest.main()
