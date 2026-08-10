#!/usr/bin/env python3
"""Build a Solv-like object directly from RDKit + AMSOLcpp results.

This script avoids writing an intermediate .solv file. It converts a 3D RDKit
molecule into AMSOLcpp atom tuples, runs both water and hexadecane calculations,
and constructs a Solv-compatible object with the same public attributes used by
traditional .solv parsers.
"""

from __future__ import annotations

import argparse
import sys
from typing import Any, Iterable, List, Optional, Sequence, Tuple

import amsolcpp
from rdkit import Chem
from rdkit.Chem import AllChem


class SolvError(RuntimeError):
    """Raised when the in-memory Solv object cannot be assembled."""


class MultiSolvException(ValueError):
    """Raised for malformed multi-record solv data."""


class Solv(object):
    """Reads .solv files from AMSOL output, or builds from in-memory rows."""

    def __init__(self, solvFileName: Optional[str] = None, rows: Optional[Sequence[Sequence[float]]] = None):
        self.name = "fake"
        self.charge: List[float] = []
        self.polarSolv: List[float] = []
        self.apolarSolv: List[float] = []
        self.solv: List[float] = []
        self.surface: List[float] = []
        self.totalAtoms: int = 0
        self.totalCharge: float = 0.0
        self.totalPolarSolv: float = 0.0
        self.totalSurface: float = 0.0
        self.totalApolarSolv: float = 0.0
        self.totalSolv: float = 0.0

        if solvFileName is not None:
            self._read_solv_file(solvFileName)
        elif rows is not None:
            self._load_rows(rows)

    def _read_solv_file(self, solvFileName: str) -> None:
        with open(solvFileName, "r", encoding="utf-8") as solvfile:
            try:
                for line in solvfile:
                    tokens = line.split()
                    if self.name == "fake":
                        self.name = tokens[0]
                        self.totalAtoms = int(tokens[1])
                        self.totalCharge = float(tokens[2])
                        self.totalPolarSolv = float(tokens[3])
                        self.totalSurface = float(tokens[4])
                        self.totalApolarSolv = float(tokens[5])
                        self.totalSolv = float(tokens[6])
                    else:
                        try:
                            self.charge.append(float(tokens[0]))
                            self.polarSolv.append(float(tokens[1]))
                            self.surface.append(float(tokens[2]))
                            self.apolarSolv.append(float(tokens[3]))
                            self.solv.append(float(tokens[4]))
                        except ValueError as exc:
                            raise MultiSolvException(line) from exc
            except (StopIteration, MultiSolvException):
                pass

    def _load_rows(self, rows: Sequence[Sequence[float]]) -> None:
        if not rows:
            raise SolvError("rows must not be empty")
        if len(rows) < 1:
            raise SolvError("expected at least one atom row")
        for row in rows:
            if len(row) != 5:
                raise SolvError("each row must contain five values: charge, polar, surface, apolar, solv")
            self.charge.append(float(row[0]))
            self.polarSolv.append(float(row[1]))
            self.surface.append(float(row[2]))
            self.apolarSolv.append(float(row[3]))
            self.solv.append(float(row[4]))
        self.totalAtoms = len(self.charge)
        self.totalCharge = 0.0
        self.totalPolarSolv = sum(self.polarSolv)
        self.totalSurface = sum(self.surface)
        self.totalApolarSolv = sum(self.apolarSolv)
        self.totalSolv = sum(self.solv)



def mol_to_amsol_atoms(mol: Chem.Mol) -> List[Tuple[int, float, float, float]]:
    if mol.GetNumConformers() == 0:
        mol = Chem.AddHs(mol)
        try:
            AllChem.EmbedMolecule(mol, randomSeed=0xF00D)
        except Exception as exc:
            raise SolvError(f"failed to generate a 3D conformer: {exc}") from exc

    if mol.GetNumConformers() == 0:
        raise SolvError("molecule has no conformers after embedding")

    conf = mol.GetConformer()
    atoms: List[Tuple[int, float, float, float]] = []
    for atom in mol.GetAtoms():
        pos = conf.GetAtomPosition(atom.GetIdx())
        atoms.append((atom.GetAtomicNum(), pos.x, pos.y, pos.z))
    return atoms


def build_solv_from_rdkit(
    mol: Chem.Mol,
    name: Optional[str] = None,
    charge: Optional[int] = None,
    verbose: bool = False,
) -> Solv:
    """Build a Solv object from an RDKit Mol with 3-D coordinates.

    The RDKit mol is passed directly to the C++ binding which handles both
    atom extraction (via ToBinary/MolPickler) and the full dual-solvent
    AM1/CM2/SM5.42R calculation in one call.  Per-field arrays returned by
    the C++ SolvDescriptors object are assigned straight to Solv attributes
    — no list comprehensions, no intermediate Python atom loop.

    The molecule must already have explicit H atoms and at least one 3-D
    conformer.  If no conformer is present, a SolvError is raised.
    """
    if name is None:
        name = mol.GetProp("_Name") if mol.HasProp("_Name") else "mol"

    if charge is None:
        charge = Chem.GetFormalCharge(mol)

    options = amsolcpp.CalculationOptions()
    options.output_decimal_precision = 2
    options.molecular_charge = charge

    try:
        desc = amsolcpp.calculate_solv_descriptors_from_rdkit_mol(mol, options)
    except Exception as exc:
        raise SolvError(f"AMSOL dual-solvent calculation failed: {exc}") from exc

    if not desc.converged_water or not desc.converged_hexadecane:
        raise SolvError(
            f"AMSOL calculation did not converge: "
            f"water={desc.converged_water}, hex={desc.converged_hexadecane}"
        )

    if verbose:
        print("cs_coeff:", desc.cs_coeff)

    # Populate Solv directly from the C++ per-field arrays — no comprehensions.
    solv = Solv()
    solv.name = name
    solv.charge = list(desc.charges)
    solv.polarSolv = list(desc.polar_diffs)
    solv.surface = list(desc.surfaces)
    solv.apolarSolv = list(desc.apolar_diffs)
    solv.solv = list(desc.solv_diffs)
    solv.totalAtoms = len(solv.charge)
    solv.totalCharge = float(charge)
    solv.totalPolarSolv = desc.total_diff_polar
    solv.totalSurface = desc.total_surface
    solv.totalApolarSolv = desc.total_diff_apolar
    solv.totalSolv = desc.total_solv_diff
    return solv


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description="Build a Solv-like object directly from RDKit and AMSOLcpp")
    parser.add_argument("smiles", help="SMILES string or RDKit mol block")
    parser.add_argument("--name", default=None, help="Name to attach to the Solv object")
    parser.add_argument("--charge", type=int, default=None, help="Formal molecular charge")
    parser.add_argument("--verbose", action="store_true", help="Print intermediate diagnostics")
    args = parser.parse_args(argv)

    try:
        mol = Chem.MolFromSmiles(args.smiles)
        if mol is None:
            raise SolvError("RDKit could not parse the SMILES string")
        mol = Chem.AddHs(mol)
        solv = build_solv_from_rdkit(mol, name=args.name, charge=args.charge, verbose=args.verbose)
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    print(f"name={solv.name}")
    print(f"total_atoms={solv.totalAtoms}")
    print(f"total_charge={solv.totalCharge}")
    print(f"total_polar={solv.totalPolarSolv}")
    print(f"total_surface={solv.totalSurface}")
    print(f"total_apolar={solv.totalApolarSolv}")
    print(f"total_solv={solv.totalSolv}")
    print("atom_rows=")
    for idx in range(solv.totalAtoms):
        print(idx, solv.charge[idx], solv.polarSolv[idx], solv.surface[idx], solv.apolarSolv[idx], solv.solv[idx])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
