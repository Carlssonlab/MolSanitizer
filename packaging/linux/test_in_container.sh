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
#   MSANI_TEST_RDKIT exact runtime RDKit version     (optional)
set -euo pipefail

tag=${MSANI_TAG:?Set MSANI_TAG to a CPython tag such as cp312}
wheels=${MSANI_WHEELS:-/wheels}
tests=${MSANI_TESTS:-/test-input}
rdkit=${MSANI_TEST_RDKIT:-}

# Resolve the wheel before changing directories. Never fall back to a source
# installation or accept multiple wheels when checking one artifact.
wheels=$(cd "$wheels" && pwd)
shopt -s nullglob
wheel_paths=("$wheels"/*"$tag"*.whl)
if [ "${#wheel_paths[@]}" -ne 1 ]; then
    echo "Expected exactly one $tag wheel in $wheels" >&2
    exit 1
fi

work=$(mktemp -d /tmp/msani-test-XXXXXXXX)
trap 'rm -rf "$work"' EXIT
cp -a "$tests/." "$work/"
cd /tmp

# Apply the same constraint to core and optional dependencies, so an extras
# installation cannot silently upgrade the selected runtime RDKit version.
constraints=()
if [ -n "$rdkit" ]; then
    printf 'rdkit==%s\n' "$rdkit" > "$work/constraints.txt"
    constraints=(--constraint "$work/constraints.txt")
fi

python -m pip install --only-binary=rdkit "${constraints[@]}" "${wheel_paths[0]}"
python -m pip check

# Verify the core installation before adding the optional PDBQT dependency.
python - <<'PY'
from pathlib import Path
import rdkit
from rdkit import rdBase
import msani, msani_confgen_cpp, msani_stereoisomers, amsolcpp
from amsolcpp import _amsolcpp
for module in (rdkit, msani, msani_confgen_cpp, msani_stereoisomers, amsolcpp, _amsolcpp):
    path = Path(module.__file__).resolve()
    assert 'site-packages' in path.parts, path
    print(module.__name__, path)
print('Runtime RDKit:', rdBase.rdkitVersion)
PY

python -m pip install --only-binary=rdkit "${constraints[@]}" "${wheel_paths[0]}[pdbqt]"
python -m pip check
python -m pip freeze

python - <<'PY'
from importlib.metadata import version
import os
import re
from rdkit import rdBase

expected = os.environ.get('MSANI_TEST_RDKIT')
installed = version('rdkit')
if expected and installed != expected:
    raise SystemExit(f'Expected Python RDKit {expected}, found {installed}')
def base_version(value):
    match = re.match(r'^(\d+)\.(\d+)\.(\d+)', value)
    if match is None:
        raise SystemExit(f'Unrecognized RDKit version: {value}')
    return tuple(map(int, match.groups()))

if base_version(rdBase.rdkitVersion) != base_version(installed):
    raise SystemExit(f'RDKit module {rdBase.rdkitVersion} disagrees with installed metadata {installed}')
print(f'Verified runtime RDKit: {installed}')
PY

python -m unittest discover -s "$work" -p 'test*.py' -v
