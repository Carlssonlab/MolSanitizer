# Third-party binary notices

MolSanitizer's Linux wheels statically link RDKit 2025.09.5. `auditwheel` may
also bundle the shared Boost.Serialization, Boost.Iostreams, bzip2, liblzma,
and zstd libraries used by the build toolchain. Their license texts are kept
in this directory and included in wheel metadata.

The exact libraries in a release remain determined by `auditwheel show` and
the repaired wheel's `molsanitizer.libs` directory. When that list changes,
update these notices before publishing.

macOS arm64 wheels use RDKit 2025.09.1 and `delocate` to bundle RDKit,
Boost.Serialization, Boost.Iostreams, bzip2, liblzma, zstd, zlib, and libc++.
The zlib and libc++ license texts are also included in wheel metadata.
Inspect `delocate-listdeps` and the repaired wheel's `.dylibs` directory for
the exact library set; `packaging/macos/check_wheel.py` verifies that native
libraries are arm64 and reference only wheel-local or macOS system libraries.
