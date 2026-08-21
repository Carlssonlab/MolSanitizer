"""Build an in-memory Solv object from an RDKit molecule using AMSOLcpp."""

from __future__ import annotations

from typing import List, Optional

import amsolcpp
from rdkit import Chem


class SolvError(RuntimeError):
    """Raised when AMSOLcpp cannot produce a usable solvation result."""


class Solv:
    """Solvation values consumed by the DB2 layout writer."""

    def __init__(self):
        self.name = "fake"
        self.charge: List[float] = []
        self.polarSolv: List[float] = []
        self.apolarSolv: List[float] = []
        self.solv: List[float] = []
        self.surface: List[float] = []
        self.totalAtoms = 0
        self.totalCharge = 0.0
        self.totalPolarSolv = 0.0
        self.totalSurface = 0.0
        self.totalApolarSolv = 0.0
        self.totalSolv = 0.0


def build_solv_from_rdkit(
    mol: Chem.Mol,
    name: Optional[str] = None,
    charge: Optional[int] = None,
    verbose: bool = False,
) -> Solv:
    """Calculate dual-solvent descriptors without intermediate files."""
    if mol.GetNumConformers() == 0:
        raise SolvError("molecule must have a 3-D conformer")

    if name is None:
        name = mol.GetProp("_Name") if mol.HasProp("_Name") else "mol"
    if charge is None:
        charge = Chem.GetFormalCharge(mol)

    options = amsolcpp.CalculationOptions()
    options.output_decimal_precision = 2
    options.molecular_charge = charge

    try:
        descriptors = amsolcpp.calculate_solv_descriptors_from_rdkit_mol(
            mol, options,
        )
    except Exception as exc:
        raise SolvError(
            f"AMSOLcpp dual-solvent calculation failed: {exc}"
        ) from exc

    if not descriptors.converged_water or not descriptors.converged_hexadecane:
        raise SolvError(
            "AMSOLcpp calculation did not converge: "
            f"water={descriptors.converged_water}, "
            f"hex={descriptors.converged_hexadecane}"
        )

    if verbose:
        print("cs_coeff:", descriptors.cs_coeff)

    solv = Solv()
    solv.name = name
    solv.charge = list(descriptors.charges)
    solv.polarSolv = list(descriptors.polar_diffs)
    solv.surface = list(descriptors.surfaces)
    solv.apolarSolv = list(descriptors.apolar_diffs)
    solv.solv = list(descriptors.solv_diffs)
    solv.totalAtoms = len(solv.charge)
    solv.totalCharge = float(charge)
    solv.totalPolarSolv = descriptors.total_diff_polar
    solv.totalSurface = descriptors.total_surface
    solv.totalApolarSolv = descriptors.total_diff_apolar
    solv.totalSolv = descriptors.total_solv_diff
    return solv
