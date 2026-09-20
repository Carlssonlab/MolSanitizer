# RDKit static SDK notice

MolSanitizer binary wheels are intended to link a private, static RDKit SDK
into their native extensions. The default SDK is RDKit **2025.09.1**
(`Release_2025_09_1`). The static archives and headers are build inputs; they
are not committed to this repository and are not copied into the wheel as
unused files.

The source incorporated by that link is covered by the BSD 3-Clause license in
[`LICENSE`](LICENSE). A distribution build must also preserve the notices for
RDKit's enabled third-party dependencies (for example Boost) from the exact
SDK build. The packaging image is the source of record for those build inputs.

Builds install `rdkit==2025.9.1` for Python-level build checks. The wheel's
runtime dependency remains `rdkit>=2025.3.1`; compatibility must be tested
separately from the SDK version. Native
code exchanges serialized molecules and Python values at that boundary; it
must not exchange RDKit C++ object pointers with the separately installed
Python RDKit package.

Use `-DMSANI_RDKIT_LINKAGE=static` and point `CMAKE_PREFIX_PATH` at the pinned
static SDK for a redistributable Linux wheel. `shared` remains the local
development default and must not be uploaded as a portable wheel.
