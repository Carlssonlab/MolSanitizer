#!/usr/bin/env bash
# Run on the Docker host, not inside a container.
#
# The build and test commands live in build_in_container.sh and
# test_in_container.sh; this wrapper only supplies containers and mounts, so
# local wheels and CI wheels come from the same code.
#
# Environment:
#   MSANI_BUILD_IMAGE     builder image           (default msani-builder:rdkit-2025.09.1)
#   MSANI_BUILD_JOBS      compile parallelism     (default 2)
#   MSANI_BUILD_RDKIT     PyPI RDKit version      (default 2025.9.1)
#   MSANI_PYTHON_VERSIONS space separated list    (default 3.10 .. 3.14)
#   MSANI_DOCKER_PLATFORM Docker platform         (default linux/amd64)
set -euo pipefail

repo=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
image=${MSANI_BUILD_IMAGE:-msani-builder:rdkit-2025.09.1}
jobs=${MSANI_BUILD_JOBS:-2}
rdkit=${MSANI_BUILD_RDKIT:-2025.9.1}
read -r -a versions <<< "${MSANI_PYTHON_VERSIONS:-3.10 3.11 3.12 3.13 3.14}"
platform=${MSANI_DOCKER_PLATFORM:-linux/amd64}
docker info >/dev/null
docker image inspect "$image" >/dev/null
mkdir -p "$repo/dist"
out=$(mktemp -d "$repo/dist/linux-wheels-XXXXXXXX")
mkdir -p "$out/wheels" "$out/logs"
printf 'Output: %s\n' "$out"
failed=0

for version in "${versions[@]}"; do
    tag="cp${version//./}"
    printf '\nBuilding Python %s\n' "$version"
    if docker run --rm --platform "$platform" \
        --mount "type=bind,source=$repo,target=/source,readonly" \
        --mount "type=bind,source=$out,target=/out" \
        -e "MSANI_TAG=$tag" -e "MSANI_BUILD_JOBS=$jobs" \
        -e "MSANI_BUILD_RDKIT=$rdkit" \
        -e "MSANI_UID=$(id -u)" -e "MSANI_GID=$(id -g)" \
        "$image" bash /source/packaging/linux/build_in_container.sh \
        2>&1 | tee "$out/logs/$tag-build.log"
    then
        printf '\nTesting Python %s in a clean container\n' "$version"
    else
        printf '%s BUILD_FAILED\n' "$version" >> "$out/status.txt"
        failed=1
        continue
    fi

    if docker run --rm --platform "$platform" \
        --mount "type=bind,source=$out/wheels,target=/wheels,readonly" \
        --mount "type=bind,source=$repo,target=/source,readonly" \
        -e "MSANI_TAG=$tag" -e "MSANI_TESTS=/source/test" \
        "python:$version-slim" bash /source/packaging/linux/test_in_container.sh \
        2>&1 | tee "$out/logs/$tag-test.log"
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
