# Local macOS wheels

Run natively on the Mac whose architecture you are targeting: Apple Silicon
produces arm64 wheels, an Intel Mac produces x86_64 wheels. Cross-architecture
builds are not supported, because every wheel is import-tested and run against
the suite before it is accepted. Needs native Conda for that architecture and
Apple's Xcode Command Line Tools (`xcode-select --install`). Network access is needed for the SDK,
standalone CPython interpreters, and Python packages.

From the repository root:

```bash
bash packaging/macos/build_wheels.sh
```

The wrapper provisions an isolated Conda C++ SDK with RDKit 2025.09.1, Boost
1.86, and Eigen 3.4 under `build/macos-<arch>/sdk`. It installs uv into a local tool
venv and uses uv-managed CPython to build and test Python 3.10–3.14. Existing
Conda environments are not changed. SDK packages are recorded in
`build/macos-<arch>/sdk-explicit.txt`; caches and interpreters stay in `build/macos-<arch>`.

For a smaller matrix or more build workers:

```bash
MSANI_BUILD_JOBS=4 bash packaging/macos/build_wheels.sh --versions 3.12
```

To use an existing arm64 SDK and uv installation directly:

```bash
python packaging/macos/build_wheels.py \
  --sdk /absolute/path/to/sdk --rdkit-version 2025.9.1 \
  --uv /absolute/path/to/uv --versions 3.12 3.13
```

The SDK needs RDKit's CMake config, libraries, headers, and matching Boost and
Eigen development packages. The PyPI RDKit version must match the SDK.
`MSANI_MACOS_SDK`, `MSANI_BUILD_RDKIT`, `MSANI_BUILD_JOBS`, and
`MSANI_BOOTSTRAP_PYTHON` override the wrapper defaults.

Each run creates `dist/macos-<arch>-wheels-*/` containing repaired wheels,
build/test logs, build settings, `status.txt`, and `SHA256SUMS`. The script
returns nonzero if any version fails. Check `status.txt` before distributing:
a repaired wheel can still be present when its installation tests fail.

Builds use fresh temporary directories and explicitly target arm64. The
[delocate](https://github.com/matthew-brett/delocate) repair step bundles
required non-system dylibs, rewrites library paths, checks arm64 compatibility,
and sanitizes RPATHs. macOS 11.0 is the requested deployment target; delocate
may raise the wheel's minimum macOS tag to match bundled SDK libraries.
Do not manually lower the final tag. Use `--deployment-target MAJOR.MINOR`
when intentionally targeting a newer macOS version.

Each repaired wheel is installed in a fresh standalone CPython environment,
with no Conda, SDK, Homebrew, or DYLD search paths. Validation imports every
native extension, checks package metadata and dependencies, installs the
`pdbqt` extra, and runs the repository unittest suite outside the source tree.
Testing on this Mac does not establish support for older macOS releases;
validate on the oldest supported OS before a release. The scripts do not
upload or publish anything.
