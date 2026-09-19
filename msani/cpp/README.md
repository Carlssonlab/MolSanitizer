# Native code map

All C++ source lives here. Python workflow code remains in `msani/conformers`,
`msani/moltransform`, and `msani/amsol`.

| Directory | What to edit here |
| --- | --- |
| `conformers/` | Embedding, RMSD, and stochastic sampling algorithms |
| `stereoisomers/` | Stereoisomer enumeration and stereo flippers |
| `amsol/include/amsolcpp/` | Vendored AMSOL public C++ interfaces |
| `amsol/src/` | Vendored AMSOL implementation |
| `bindings/` | Python argument conversion and pybind11 registrations |

Headers stay beside their implementation, except AMSOL, whose upstream
`include/` and `src/` structure is preserved for updates. Keep AMSOL's license
and upstream documentation in `amsol/`.

## Build structure

The repository root CMake file loads this directory. Each component CMake file
selects its implementation and corresponding file from `bindings/`.
Shared build helpers are in the repository's `cmake/` directory.

This source reorganization preserves three extension modules:

- `msani_confgen_cpp`: `bindings/conformers.cpp`
- `msani_stereoisomers`: `bindings/stereoisomers.cpp`
- `amsolcpp._amsolcpp`: `bindings/amsol.cpp`

The Python `amsolcpp` wrapper is maintained at `msani/amsol/amsolcpp/` and
installed as a top-level package by CMake. Public imports are unchanged.
AMSOL uses RDKit's Python API; the other two extensions link RDKit C++.

Combining these into a single `msani._native` extension is a separate binding
and API migration. The common source directory alone does not consolidate
their binary linkage.

Build from the repository root using the normal pip or CMake workflow.
Use a fresh build directory after moving source paths, or reconfigure an
existing directory before building. See the root `build_rdkit_SDK.md` for
the static SDK wheel workflow.

## Distribution

The source archive includes this tree. Binary wheels exclude `msani/cpp/`;
CMake installs the compiled modules and their `.pyi` files separately.
License notices are included through the wheel's license metadata.
