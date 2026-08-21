"""Stable public entry point for in-memory MOL2-to-DB2 conversion."""

from __future__ import annotations

from os import PathLike
from typing import Any

from msani.db2.hierarchy import build_conformer_hierarchy, serialize_db2
from msani.db2.molecule import MoleculeData, prepare_molecule_for_db2


def mol2db2(
    mol2data: MoleculeData,
    solvdata: Any,
    clashfile: str | PathLike[str] | None = None,
    disttol: float = 0.001,
) -> str:
    """Convert in-memory molecular and solvation data to DB2 text.

    ``mol2data`` is prepared in place, preserving the historical metadata and
    atom-property side effects. Coordinates and floating-point calculations
    are not reordered during hierarchy construction or serialization.
    """

    _validate_legacy_clash_file(clashfile)
    prepared_molecule = prepare_molecule_for_db2(mol2data)
    conformer_hierarchy = build_conformer_hierarchy(
        prepared_molecule, tolerance=disttol
    )
    return serialize_db2(prepared_molecule, solvdata, conformer_hierarchy)


def _validate_legacy_clash_file(
    clashfile: str | PathLike[str] | None,
) -> None:
    """Preserve validation of the legacy clash-file option.

    Clash filtering was disabled in the pre-refactor hierarchy builder, but a
    supplied file was still opened and parsed. Keeping that input boundary
    avoids changing errors observed by callers without suggesting that these
    rules affect generated conformer sets.
    """

    if clashfile is None:
        return
    with open(clashfile) as clash_rules:
        for line in clash_rules:
            tokens = line.split()
            _constraint = tokens[0]
            int(tokens[1])
            int(tokens[2])
            _first_atom_type = tokens[3]
            _second_atom_type = tokens[4]
            float(tokens[5])
