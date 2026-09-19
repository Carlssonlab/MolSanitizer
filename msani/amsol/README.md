# Python solvation integration

`rdkit_amsol_to_solv.py` adapts AMSOL calculation results for MolSanitizer's
DB2 output. `amsolcpp/__init__.py` provides the Python wrapper installed as
the top-level `amsolcpp` package.

The native source, license, and scientific documentation are in
[`../cpp/amsol/`](../cpp/amsol/README.md). Python bindings are in
[`../cpp/bindings/amsol.cpp`](../cpp/bindings/amsol.cpp).
