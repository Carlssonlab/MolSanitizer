#!/usr/bin/env bash
# Run natively on an Apple Silicon Mac with Conda and Xcode Command Line Tools.
set -euo pipefail
repo=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
if [[ $(uname -s) != Darwin || $(uname -m) != arm64 ]]; then
    printf 'Run on native macOS arm64, outside Rosetta.\n' >&2
    exit 1
fi
xcrun --find clang++ >/dev/null
conda_exe=${CONDA_EXE:-$(command -v conda || true)}
if [[ -z "$conda_exe" || ! -x "$conda_exe" ]]; then
    printf 'Install native arm64 Conda, or set CONDA_EXE to its executable.\n' >&2
    exit 1
fi
bootstrap_python=${MSANI_BOOTSTRAP_PYTHON:-$("$conda_exe" info --base)/bin/python}
cache="$repo/build/macos"
sdk=${MSANI_MACOS_SDK:-$cache/sdk}
rdkit=${MSANI_BUILD_RDKIT:-2025.9.1}
mkdir -p "$cache"
# Keep all new environments and caches local to the repository.
export CONDA_PKGS_DIRS="$cache/conda-pkgs"
export PIP_CACHE_DIR="$cache/pip-cache"
if [[ ! -f "$sdk/lib/cmake/rdkit/rdkit-config.cmake" ]]; then
    "$conda_exe" create --yes --prefix "$sdk" --override-channels \
        -c conda-forge --platform osx-arm64 "librdkit-dev=$rdkit" \
        'libboost-devel=1.86' 'eigen=3.4'
fi
"$conda_exe" list --prefix "$sdk" --explicit > "$cache/sdk-explicit.txt"
if [[ ! -x "$cache/tools/bin/uv" ]]; then
    "$bootstrap_python" -m venv "$cache/tools"
    "$cache/tools/bin/python" -m pip install uv
fi
exec "$cache/tools/bin/python" "$repo/packaging/macos/build_wheels.py" \
    --uv "$cache/tools/bin/uv" --sdk "$sdk" --rdkit-version "$rdkit" \
    --jobs "${MSANI_BUILD_JOBS:-2}" "$@"
