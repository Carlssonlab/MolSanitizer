#!/usr/bin/env bash
# Provision the conda-forge RDKit C++ SDK for this Mac's architecture.
#
# Single source of truth for the macOS SDK: build_wheels.sh calls this before
# building, and the GitHub Actions SDK job calls it before packing the prefix
# into a release asset. Prints the SDK prefix on the last line.
#
# Environment:
#   MSANI_MACOS_SDK   SDK prefix override
#   MSANI_BUILD_RDKIT PyPI/conda RDKit version (default 2025.9.1)
#   CONDA_EXE         conda executable
set -euo pipefail

repo=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
if [[ $(uname -s) != Darwin ]]; then
    printf 'Run on macOS.\n' >&2
    exit 1
fi
arch=$(uname -m)
case "$arch" in
    arm64) subdir=osx-arm64 ;;
    x86_64) subdir=osx-64 ;;
    *) printf 'Unsupported macOS architecture: %s\n' "$arch" >&2; exit 1 ;;
esac
conda_exe=${CONDA_EXE:-$(command -v conda || true)}
if [[ -z "$conda_exe" || ! -x "$conda_exe" ]]; then
    printf 'Install native %s Conda, or set CONDA_EXE to its executable.\n' "$arch" >&2
    exit 1
fi
cache="$repo/build/macos-$arch"
sdk=${MSANI_MACOS_SDK:-$cache/sdk}
rdkit=${MSANI_BUILD_RDKIT:-2025.9.1}
mkdir -p "$cache"
# Keep all new environments and caches local to the repository.
export CONDA_PKGS_DIRS="$cache/conda-pkgs"

# The SDK supplies RDKit's C++ libraries, headers and CMake exports only; the
# Python RDKit API is installed separately by the wheel build.
if [[ ! -f "$sdk/lib/cmake/rdkit/rdkit-config.cmake" ]]; then
    "$conda_exe" create --yes --prefix "$sdk" --override-channels \
        -c conda-forge --platform "$subdir" "librdkit-dev=$rdkit" \
        'libboost-devel=1.86' 'eigen=3.4'
fi
"$conda_exe" list --prefix "$sdk" --explicit > "$cache/sdk-explicit.txt"
printf '%s\n' "$sdk"
