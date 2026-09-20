#!/usr/bin/env bash
# Build and auditwheel-repair one CPython wheel inside a manylinux_2_28
# environment that already carries the RDKit SDK.
#
# Used both by build_wheels.sh (via docker run) and by the GitHub Actions
# container job, so the local and CI wheels come from identical commands.
#
# Inputs (environment):
#   MSANI_TAG           cp310 .. cp314          (required)
#   MSANI_SOURCE        repository checkout     (default /source)
#   MSANI_OUT           repaired wheel output   (default /out/wheels)
#   MSANI_BUILD_RDKIT   PyPI RDKit version matching the SDK
#   MSANI_BUILD_JOBS    compile parallelism
#   MSANI_UID/MSANI_GID chown the output back to the Docker host user
set -euo pipefail

tag=${MSANI_TAG:?Set MSANI_TAG to a CPython tag such as cp312}
source_dir=${MSANI_SOURCE:-/source}
out=${MSANI_OUT:-/out/wheels}
rdkit=${MSANI_BUILD_RDKIT:-2025.9.1}
jobs=${MSANI_BUILD_JOBS:-2}
rdkit_prefix=${RDKIT_PREFIX:-/opt/msani-rdkit}
boost_prefix=${BOOST_PREFIX:-/opt/boost}

# manylinux_2_28 is the packaging target on every supported architecture.
plat=${MSANI_PLAT:-manylinux_2_28_$(uname -m)}

mkdir -p "$out"
if [[ -n "${MSANI_UID:-}" && -n "${MSANI_GID:-}" ]]; then
    trap 'chown -R "$MSANI_UID:$MSANI_GID" "$out"' EXIT
fi

export PATH="/opt/python/$tag-$tag/bin:$PATH"
export RDBASE="$rdkit_prefix"
export CMAKE_BUILD_PARALLEL_LEVEL="$jobs"

# Use a writable copy; never reuse host CMake artifacts or modify the checkout.
work=$(mktemp -d /tmp/msani-XXXXXXXX)
mkdir -p "$work/source" "$work/raw"
tar -C "$source_dir" --exclude='./.git' --exclude='./build' --exclude='./dist' \
    --exclude='__pycache__' -cf - . | tar -C "$work/source" -xf -
cd "$work/source"

python -m pip install scikit-build-core pybind11 twine "rdkit==$rdkit"
python -m pip wheel . --no-build-isolation --no-deps \
    -Cbuild-dir="$work/native-build" \
    -Ccmake.define.MSANI_RDKIT_LINKAGE=static \
    -Ccmake.define.CMAKE_PREFIX_PATH="$rdkit_prefix;$boost_prefix" \
    --wheel-dir "$work/raw"

auditwheel show "$work"/raw/molsanitizer-*.whl
auditwheel repair --plat "$plat" --wheel-dir "$out" "$work"/raw/molsanitizer-*.whl
python -m twine check --strict "$out"/*"$tag"*.whl
