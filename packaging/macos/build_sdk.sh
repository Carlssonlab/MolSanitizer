#!/usr/bin/env bash
# Provision the conda-forge RDKit C++ SDK for this Mac's architecture.
#
# Single source of truth for the macOS SDK: build_wheels.sh calls this before
# building, and the GitHub Actions SDK job calls it before packing the prefix
# into a release asset. Prints the SDK prefix on the last line.
#
# Environment:
#   MSANI_MACOS_SDK   SDK prefix override
#   MSANI_BUILD_RDKIT PyPI/conda RDKit version (default 2026.9.2)
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
rdkit=${MSANI_BUILD_RDKIT:-2026.9.2}
mkdir -p "$cache"
# Keep all new environments and caches local to the repository.
export CONDA_PKGS_DIRS="$cache/conda-pkgs"

# The SDK supplies RDKit's C++ libraries, headers and CMake exports only; the
# Python RDKit API is installed separately by the wheel build.
# Boost and Eigen follow what conda-forge built this librdkit release against.
if [[ ! -f "$sdk/lib/cmake/rdkit/rdkit-config.cmake" ]]; then
    "$conda_exe" create --yes --prefix "$sdk" --override-channels \
        -c conda-forge --platform "$subdir" "librdkit-dev=$rdkit" \
        'libboost-devel=1.90' 'eigen=5.0'
fi
# A prefix provisioned for another RDKit pin must not be reused silently.
# Compare with leading zeros dropped: conda spells 2026.9.2 as 2026.09.2.
unpad() { sed -E 's/(^|\.)0+([0-9])/\1\2/g' <<<"$1"; }
installed=$(sed -nE 's/^set\(PACKAGE_VERSION "([0-9.]+)"\)$/\1/p' \
    "$sdk/lib/cmake/rdkit/rdkit-config-version.cmake" | head -n 1)
if [[ "$(unpad "$installed")" != "$(unpad "$rdkit")" ]]; then
    printf 'SDK at %s holds RDKit %s, not %s. Remove it or set MSANI_MACOS_SDK.\n' \
        "$sdk" "${installed:-unknown}" "$rdkit" >&2
    exit 1
fi
"$conda_exe" list --prefix "$sdk" --explicit > "$cache/sdk-explicit.txt"
printf '%s\n' "$sdk"
