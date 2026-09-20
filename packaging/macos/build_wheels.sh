#!/usr/bin/env bash
# Run natively on a Mac with Conda and Xcode Command Line Tools.
# Apple Silicon produces arm64 wheels; an Intel Mac (or a macos-15-intel
# runner) produces x86_64 wheels. Cross-architecture builds are not supported:
# every wheel is import-tested and run against the suite before acceptance.
set -euo pipefail
repo=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
if [[ $(uname -s) != Darwin ]]; then
    printf 'Run on macOS.\n' >&2
    exit 1
fi
arch=$(uname -m)
xcrun --find clang++ >/dev/null
conda_exe=${CONDA_EXE:-$(command -v conda || true)}
if [[ -z "$conda_exe" || ! -x "$conda_exe" ]]; then
    printf 'Install native %s Conda, or set CONDA_EXE to its executable.\n' "$arch" >&2
    exit 1
fi
bootstrap_python=${MSANI_BOOTSTRAP_PYTHON:-$("$conda_exe" info --base)/bin/python}
cache="$repo/build/macos-$arch"
rdkit=${MSANI_BUILD_RDKIT:-2025.9.1}
export PIP_CACHE_DIR="$cache/pip-cache"
# Provision (or reuse) the SDK, then read the prefix it reports.
sdk=$(bash "$repo/packaging/macos/build_sdk.sh" | tail -n 1)
if [[ ! -x "$cache/tools/bin/uv" ]]; then
    "$bootstrap_python" -m venv "$cache/tools"
    "$cache/tools/bin/python" -m pip install uv
fi
exec "$cache/tools/bin/python" "$repo/packaging/macos/build_wheels.py" \
    --uv "$cache/tools/bin/uv" --sdk "$sdk" --rdkit-version "$rdkit" \
    --arch "$arch" --jobs "${MSANI_BUILD_JOBS:-2}" "$@"
