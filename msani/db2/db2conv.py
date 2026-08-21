"""Public entry point for in-memory DB2 conversion."""

from __future__ import annotations

from typing import Any

from msani.db2.layout import build_db2_layout, serialize_db2
from msani.db2.molecule import MoleculeData, prepare_molecule_for_db2


def db2converter(
    molecule: MoleculeData,
    solvation_data: Any,
    disttol: float = 0.001,
) -> str:
    """Convert in-memory molecular and solvation data to DB2 text.

    ``molecule`` is prepared in place. Coordinates and floating-point
    calculations are not reordered during DB2 layout construction or
    serialization.
    """

    prepared_molecule = prepare_molecule_for_db2(molecule)
    layout = build_db2_layout(prepared_molecule, tolerance=disttol)
    return serialize_db2(prepared_molecule, solvation_data, layout)
