"""Ergonomic Python API backed entirely by the native AMSOLcpp extension."""

import math

from ._amsolcpp import (
    __version__,
    AmsolError,
    Atom,
    AtomicSolvationResult,
    CalculationOptions,
    CalculationResult,
    DiagonalizationError,
    DualSolventResult,
    InputError,
    NumericalError,
    ScfConvergenceError,
    ScfIterationRecord,
    Solvent,
    UnsupportedCalculationError,
    UnsupportedElementError,
    atoms_from_rdkit_mol,
    calculate as _native_calculate,
    calculate_batch as _native_calculate_batch,
    calculate_from_input,
    calculate_from_rdkit_mol,
    calculate_solv_descriptors_from_rdkit_mol,
    calculate_water_and_hexadecane as _native_calculate_dual,
    run_legacy_input,
)
from . import _amsolcpp as _native_module

# True when this build was compiled with RDKit support (AMSOLCPP_WITH_RDKIT).
# False means the three rdkit functions exist but raise RuntimeError on call.
rdkit_support: bool = getattr(_native_module, "rdkit_support", False)


_ATOMIC_NUMBERS = {
    "H": 1, "C": 6, "N": 7, "O": 8, "F": 9, "SI": 14,
    "P": 15, "S": 16, "CL": 17, "BR": 35, "I": 53,
}
_OPTION_FIELDS = (
    "solvent",
    "molecular_charge",
    "multiplicity",
    "max_scf_iterations",
    "legacy_scf_iterations",
    "energy_tolerance",
    "density_tolerance",
    "commutator_tolerance",
    "enable_scf_rescue",
    "collect_iteration_diagnostics",
    "output_decimal_precision",
)


def _atomic_number(value):
    if isinstance(value, bool):
        raise TypeError("atomic number must be an integer or supported symbol")
    if isinstance(value, str):
        try:
            return _ATOMIC_NUMBERS[value.strip().upper()]
        except KeyError as exc:
            raise ValueError("unsupported element symbol: {0!r}".format(value)) from exc
    if isinstance(value, int):
        converted = value
    else:
        try:
            converted = int(value)
        except (TypeError, ValueError) as exc:
            raise TypeError(
                "atomic number must be an integer or supported symbol"
            ) from exc
    if converted < -2147483648 or converted > 2147483647:
        raise OverflowError("atomic number exceeds the C++ int range")
    if converted != value:
        raise TypeError("atomic number must be an integer or supported symbol")
    return converted


def _coerce_atom(value, index):
    if isinstance(value, Atom):
        coordinates = (value.x_angstrom, value.y_angstrom, value.z_angstrom)
        number = value.atomic_number
    else:
        try:
            if len(value) != 4:
                raise ValueError
            number, *coordinates = value
        except (TypeError, ValueError) as exc:
            raise TypeError(
                "atom {0} must be Atom or (symbol_or_number, x, y, z)".format(index)
            ) from exc
    number = _atomic_number(number)
    try:
        x, y, z = (float(component) for component in coordinates)
    except (TypeError, ValueError) as exc:
        raise TypeError("atom {0} coordinates must be real numbers".format(index)) from exc
    if not all(math.isfinite(component) for component in (x, y, z)):
        raise ValueError("atom {0} coordinates must be finite".format(index))
    return Atom(number, x, y, z)


def _coerce_atoms(atoms):
    try:
        values = list(atoms)
    except TypeError as exc:
        raise TypeError("atoms must be an iterable") from exc
    return [_coerce_atom(value, index) for index, value in enumerate(values)]


def _parse_solvent(value):
    if isinstance(value, Solvent):
        return value
    normalized = str(value).strip().lower().replace("-", "_")
    if normalized == "water":
        return Solvent.WATER
    if normalized in {
        "hexadecane",
        "genorg_hexadecane_exact_audited",
        "genorg_hexadecane",
    }:
        return Solvent.HEXADECANE
    raise ValueError("solvent must be 'water' or 'hexadecane'")


def _configured_options(
    options=None,
    *,
    charge=None,
    solvent=None,
    multiplicity=None,
    max_scf_iterations=None,
    collect_iteration_diagnostics=None
):
    if options is not None and not isinstance(options, CalculationOptions):
        raise TypeError("options must be a CalculationOptions instance")
    configured = CalculationOptions()
    if options is not None:
        for field in _OPTION_FIELDS:
            setattr(configured, field, getattr(options, field))
    if charge is not None:
        if isinstance(charge, bool) or not isinstance(charge, int):
            raise TypeError("charge must be an integer")
        configured.molecular_charge = charge
    if solvent is not None:
        configured.solvent = _parse_solvent(solvent)
    if multiplicity is not None:
        if isinstance(multiplicity, bool) or not isinstance(multiplicity, int):
            raise TypeError("multiplicity must be an integer")
        configured.multiplicity = multiplicity
    if max_scf_iterations is not None:
        if isinstance(max_scf_iterations, bool) or not isinstance(max_scf_iterations, int):
            raise TypeError("max_scf_iterations must be an integer")
        if max_scf_iterations <= 0:
            raise ValueError("max_scf_iterations must be positive")
        configured.max_scf_iterations = max_scf_iterations
    if collect_iteration_diagnostics is not None:
        configured.collect_iteration_diagnostics = bool(collect_iteration_diagnostics)
    return configured


def calculate(
    atoms,
    options=None,
    *,
    charge=None,
    solvent=None,
    multiplicity=None,
    max_scf_iterations=None,
    collect_iteration_diagnostics=None
):
    """Calculate AM1/CM2/SM5.42R properties from Atom objects or 4-tuples."""
    configured = _configured_options(
        options,
        charge=charge,
        solvent=solvent,
        multiplicity=multiplicity,
        max_scf_iterations=max_scf_iterations,
        collect_iteration_diagnostics=collect_iteration_diagnostics,
    )
    return _native_calculate(_coerce_atoms(atoms), configured)


def calculate_arrays(atomic_numbers, coordinates, **kwargs):
    """Calculate from atomic numbers and a NumPy coordinates array of shape (N, 3)."""
    try:
        import numpy as np
    except ImportError as exc:
        raise ImportError("NumPy is required only for calculate_arrays") from exc
    numbers = np.asarray(atomic_numbers)
    xyz = np.asarray(coordinates)
    if numbers.ndim != 1:
        raise ValueError("atomic_numbers must have shape (N,)")
    if numbers.dtype.kind not in "iu":
        raise TypeError("atomic_numbers must have an integer dtype")
    if xyz.ndim != 2 or xyz.shape != (numbers.shape[0], 3):
        raise ValueError("coordinates must have shape (N, 3)")
    if xyz.dtype.kind != "f":
        raise TypeError("coordinates must have a floating-point dtype")
    if not np.isfinite(xyz).all():
        raise ValueError("coordinates must contain only finite values")
    atoms = [
        (int(numbers[index]), float(xyz[index, 0]), float(xyz[index, 1]), float(xyz[index, 2]))
        for index in range(numbers.shape[0])
    ]
    return calculate(atoms, **kwargs)


def calculate_water_and_hexadecane(atoms, options=None, *, charge=None, multiplicity=None):
    """Return (water_result, hexadecane_result) with shared AM1 setup."""
    configured = _configured_options(
        options, charge=charge, multiplicity=multiplicity
    )
    result = _native_calculate_dual(_coerce_atoms(atoms), configured)
    return result.water, result.hexadecane


def calculate_batch(
    molecules,
    options=None,
    *,
    charge=None,
    solvent=None,
    continue_on_error=False
):
    """Calculate molecules sequentially, preserving input order and indexed failures."""
    configured = _configured_options(options, charge=charge, solvent=solvent)
    prepared = []
    preparation_failures = {}
    for index, molecule in enumerate(molecules):
        try:
            prepared.append(_coerce_atoms(molecule))
        except Exception as exc:
            message = "batch molecule {0} failed: {1}".format(index, exc)
            if not continue_on_error:
                raise InputError(message) from exc
            prepared.append([])
            preparation_failures[index] = message
    results = _native_calculate_batch(
        prepared,
        configured,
        continue_on_error=continue_on_error or bool(preparation_failures),
    )
    for index, message in preparation_failures.items():
        failed = CalculationResult()
        failed.warnings = [message]
        results[index] = failed
    return results


def run_input(input_text):
    """Run legacy AMSOL input text and return the supported AMSOL-style block."""
    return run_legacy_input(input_text)


run_amsol_text = run_input


def run_amsol(path):
    """Run a legacy AMSOL input file and return native formatted output."""
    with open(path, "r", encoding="utf-8") as stream:
        return run_input(stream.read())


run_amsol_input_file = run_amsol


def _main():
    import argparse
    parser = argparse.ArgumentParser(description="Run native AMSOLcpp from Python")
    parser.add_argument("input")
    parser.add_argument("-o", "--output")
    args = parser.parse_args()
    text = run_amsol(args.input)
    if args.output:
        with open(args.output, "w", encoding="utf-8") as stream:
            stream.write(text)
    else:
        print(text, end="")


__all__ = [
    "__version__", "AmsolError", "Atom", "AtomicSolvationResult",
    "CalculationOptions", "CalculationResult", "DiagonalizationError",
    "DualSolventResult", "InputError", "NumericalError", "ScfConvergenceError",
    "ScfIterationRecord", "Solvent", "UnsupportedCalculationError",
    "UnsupportedElementError", "atoms_from_rdkit_mol", "calculate",
    "calculate_arrays", "calculate_batch", "calculate_from_input",
    "calculate_from_rdkit_mol", "calculate_solv_descriptors_from_rdkit_mol",
    "calculate_water_and_hexadecane", "rdkit_support", "run_amsol",
    "run_amsol_input_file", "run_amsol_text", "run_input", "run_legacy_input",
]
