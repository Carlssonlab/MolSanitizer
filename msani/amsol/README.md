# AMSOLcpp in MolSanitizer

This directory contains the AMSOLcpp source code vendored by MolSanitizer.
MolSanitizer builds the native Python extension together with its other C++
modules; users do not need to install AMSOLcpp separately.

The upstream project is available at
[isra3l/AMSOLcpp](https://github.com/isra3l/AMSOLcpp).

## Purpose

AMSOLcpp is a C++20 port of the
[AMSOL 7.1](https://comp.chem.umn.edu/amsol/) FORTRAN 77 program. It implements
the audited `AM1 1SCF SM5.42R` workflow used by MolSanitizer to calculate the
partial charges and desolvation descriptors required for DB2 output.

Calculations run directly from an in-memory RDKit molecule. MolSanitizer does
not invoke an external AMSOL executable and does not create intermediate AMSOL
input, output, or `.solv` files.

## Supported scientific scope

- Closed-shell restricted Hartree-Fock (RHF) with the AM1 Hamiltonian
- CM2 atomic charges
- SM5.42R water and the audited GENORG hexadecane parameter set
- H, C, N, O, F, Si, P, S, Cl, Br, and I
- Cartesian structured input and legacy MOPAC Z-matrix input
- Deterministic sequential execution by default

## MolSanitizer integration

The adapter in
[`rdkit_amsol_to_solv.py`](rdkit_amsol_to_solv.py) converts the native
AMSOLcpp result into the `Solv` object consumed by MolSanitizer's DB2 hierarchy
writer.

For each molecule, the adapter:

1. Requires an RDKit molecule with a three-dimensional conformer.
2. Uses the molecule's formal charge unless a charge is supplied explicitly.
3. Runs the water and hexadecane calculations from the same molecular geometry.
4. Takes CM2 atomic charges and surface areas from the hexadecane calculation
   (dielectric constant 2.06).
5. Calculates polar and apolar desolvation descriptors as the water value minus
   the corresponding hexadecane value.
6. Rejects the result if either solvent calculation does not converge.

## Build and packaging

The top-level MolSanitizer `CMakeLists.txt` adds this directory as a CMake
subdirectory. Its Python package and compiled `_amsolcpp` extension are then
installed into the MolSanitizer wheel. This directory deliberately has no
separate `pyproject.toml`.

The native core is compiled with `-fno-fast-math` and `-ffp-contract=off` on
GCC and Clang-family compilers to preserve the floating-point behavior used for
scientific parity testing. CPU-specific native optimizations and OpenMP are not
enabled in distributed builds.

CMake first looks for an RDKit CMake package and then checks common Conda and
system locations for RDKit development files. If they are unavailable, the
extension can still build, but its direct RDKit functions will report that
RDKit support was not compiled in. MolSanitizer's DB2 workflow requires that
support.

Do not keep a separate editable installation of the standalone `amsolcpp`
Python package in the same environment. Both installations provide the same
top-level import name, so the standalone editable import hook can shadow the
version bundled with MolSanitizer.

## Validation and performance

During development, a benchmark on a random sample of 10,000 Enamine molecules
reported parity with the AMSOL 7.1 reference results for every evaluated
molecule. The reported median per-molecule speedup was 27.42×, with individual
speedups ranging from 19.89× to 32.26×.

These figures describe that validation data set and environment; performance
on other hardware and molecular data sets may differ.

## Citation

If results produced through AMSOLcpp are used in research, cite AMSOL 7.1:

> Hawkins, G. D.; Giesen, D.; Lynch, G.; Chambers, C.; Rossi, I.; et al.
> *AMSOL*, version 7.1. University of Minnesota, 2004.

## License

[AMSOL 7.1](https://comp.chem.umn.edu/amsol/) is distributed under the Apache
License, Version 2.0. AMSOLcpp follows the same license, and the source vendored
here is distributed under the terms in [`LICENSE`](LICENSE).
