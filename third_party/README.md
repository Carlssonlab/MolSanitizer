# Third-party binary notices

MolSanitizer's Linux wheels statically link RDKit 2025.09.5. `auditwheel` may
also bundle the shared Boost.Serialization, Boost.Iostreams, bzip2, liblzma,
and zstd libraries used by the build toolchain. Their license texts are kept
in this directory and included in wheel metadata.

The exact libraries in a release remain determined by `auditwheel show` and
the repaired wheel's `molsanitizer.libs` directory. When that list changes,
update these notices before publishing.
