#!/usr/bin/env bash
# Run on the Docker host, not inside a container.
set -euo pipefail

repo=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
image=${MSANI_BUILD_IMAGE:-msani-builder:rdkit-2025.09.1}
jobs=${MSANI_BUILD_JOBS:-2}
rdkit=${MSANI_BUILD_RDKIT:-2025.9.1}
docker info >/dev/null
docker image inspect "$image" >/dev/null
mkdir -p "$repo/dist"
out=$(mktemp -d "$repo/dist/linux-wheels-XXXXXXXX")
mkdir -p "$out/wheels" "$out/logs"
printf 'Output: %s\n' "$out"
failed=0

for version in 3.10 3.11 3.12 3.13 3.14; do
    tag="cp${version//./}"
    printf '\nBuilding Python %s\n' "$version"
    if docker run --rm -i --platform linux/amd64 \
        --mount "type=bind,source=$repo,target=/source,readonly" \
        --mount "type=bind,source=$out,target=/out" \
        -e "MSANI_TAG=$tag" -e "MSANI_BUILD_JOBS=$jobs" \
        -e "MSANI_BUILD_RDKIT=$rdkit" \
        -e "MSANI_UID=$(id -u)" -e "MSANI_GID=$(id -g)" \
        "$image" bash -s 2>&1 <<'BUILD' | tee "$out/logs/$tag-build.log"
set -euo pipefail
trap 'chown -R "$MSANI_UID:$MSANI_GID" /out/wheels' EXIT
export PATH="/opt/python/$MSANI_TAG-$MSANI_TAG/bin:$PATH"
export RDBASE=/opt/msani-rdkit
export CMAKE_BUILD_PARALLEL_LEVEL="$MSANI_BUILD_JOBS"
# Use a writable copy; never reuse host CMake artifacts or modify the checkout.
mkdir -p /tmp/source /tmp/raw
tar -C /source --exclude='./.git' --exclude='./build' --exclude='./dist' \
    --exclude='__pycache__' -cf - . | tar -C /tmp/source -xf -
cd /tmp/source
python -m pip install scikit-build-core pybind11 twine "rdkit==$MSANI_BUILD_RDKIT"
python -m pip wheel . --no-build-isolation --no-deps \
    -Cbuild-dir=/tmp/native-build \
    -Ccmake.define.MSANI_RDKIT_LINKAGE=static \
    -Ccmake.define.CMAKE_PREFIX_PATH='/opt/msani-rdkit;/opt/boost' \
    --wheel-dir /tmp/raw
auditwheel show /tmp/raw/molsanitizer-*.whl
auditwheel repair --plat manylinux_2_28_x86_64 \
    --wheel-dir /out/wheels /tmp/raw/molsanitizer-*.whl
python -m twine check --strict /out/wheels/*"$MSANI_TAG"*.whl
BUILD
    then
        printf '\nTesting Python %s in a clean container\n' "$version"
    else
        printf '%s BUILD_FAILED\n' "$version" >> "$out/status.txt"
        failed=1
        continue
    fi

    if docker run --rm -i --platform linux/amd64 \
        --mount "type=bind,source=$out/wheels,target=/wheels,readonly" \
        --mount "type=bind,source=$repo/test,target=/test-input,readonly" \
        -e "MSANI_TAG=$tag" \
        "python:$version-slim" bash -s 2>&1 <<'TEST' | tee "$out/logs/$tag-test.log"
set -euo pipefail
mkdir -p /tmp/tests
cp -a /test-input/. /tmp/tests/
cd /tmp
python -m pip install /wheels/*"$MSANI_TAG"*.whl
python -m pip check
# Verify the core installation before adding the optional PDBQT test dependency.
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
python -m unittest discover -s /tmp/tests -p 'test*.py' -v
TEST
    then
        printf '%s PASSED\n' "$version" >> "$out/status.txt"
    else
        printf '%s TEST_FAILED\n' "$version" >> "$out/status.txt"
        failed=1
    fi
done

(cd "$out/wheels" && sha256sum ./*.whl) > "$out/SHA256SUMS" || failed=1
printf '\nResults: %s\n' "$out"
cat "$out/status.txt"
exit "$failed"
