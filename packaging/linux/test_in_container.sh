#!/usr/bin/env bash
# Install one repaired wheel into a clean CPython image and run the suite.
#
# Deliberately run in a bare python:X.Y-slim environment, not the build image:
# a wheel that needs manylinux build-image libraries is not redistributable.
#
# Inputs (environment):
#   MSANI_TAG      cp310 .. cp314                    (required)
#   MSANI_WHEELS   directory holding repaired wheels (default /wheels)
#   MSANI_TESTS    repository test directory         (default /test-input)
set -euo pipefail

tag=${MSANI_TAG:?Set MSANI_TAG to a CPython tag such as cp312}
wheels=${MSANI_WHEELS:-/wheels}
tests=${MSANI_TESTS:-/test-input}

work=$(mktemp -d /tmp/msani-test-XXXXXXXX)
cp -a "$tests/." "$work/"
cd /tmp

python -m pip install "$wheels"/*"$tag"*.whl
python -m pip check

# Verify the core installation before adding the optional PDBQT dependency.
python - <<'PY'
from pathlib import Path
from rdkit import rdBase
import msani, msani_confgen_cpp, msani_stereoisomers, amsolcpp
from amsolcpp import _amsolcpp
for module in (msani, msani_confgen_cpp, msani_stereoisomers, amsolcpp, _amsolcpp):
    path = Path(module.__file__).resolve()
    assert 'site-packages' in path.parts, path
    print(module.__name__, path)
print('Runtime RDKit:', rdBase.rdkitVersion)
PY

python -m pip install 'meeko>=0.7.1' 'scipy'
python -m pip check
python -m pip freeze
python -m unittest discover -s "$work" -p 'test*.py' -v
