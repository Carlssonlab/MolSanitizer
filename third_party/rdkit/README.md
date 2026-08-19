# RDKit static SDK notice

MolSanitizer binary wheels are intended to link a private, static RDKit SDK
into their native extensions. The initial supported SDK is RDKit **2025.09.5**
(`Release_2025_09_5`). The static archives and headers are build inputs; they
are not committed to this repository and are not copied into the wheel as
unused files.

The source incorporated by that link is covered by the BSD 3-Clause license in
[`LICENSE`](LICENSE). A distribution build must also preserve the notices for
RDKit's enabled third-party dependencies (for example Boost) from the exact
SDK build. The packaging image is the source of record for those build inputs.

The wheel continues to install the `rdkit==2025.9.5` Python package. Native
code exchanges serialized molecules and Python values at that boundary; it
must not exchange RDKit C++ object pointers with the separately installed
Python RDKit package.

Use `-DMSANI_RDKIT_LINKAGE=static` and point `CMAKE_PREFIX_PATH` at the pinned
static SDK for a redistributable Linux wheel. `shared` remains the local
development default and must not be uploaded as a portable wheel.
