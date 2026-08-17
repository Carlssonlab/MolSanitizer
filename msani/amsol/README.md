# AMSOLcpp adapter

MolSanitizer calculates DB2 partial charges and solvation descriptors through
the [AMSOLcpp](https://github.com/isra3l/AMSOLcpp) Python binding.
`rdkit_amsol_to_solv.py` adapts the in-memory result to the object consumed by
the DB2 hierarchy writer.

The legacy AMSOL 7.1 executables, Z-matrix converter, and temporary `.solv`
workflow have been removed.
