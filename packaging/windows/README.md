# Windows wheels from a conda Command Prompt

The default workflow now uses **conda**, not the `py` launcher. It builds
CPython 3.10–3.14 x64 wheels locally with VS Build Tools. No Docker or GitHub
Actions is required.

Use `build_sdk.cmd` and `build_wheels.cmd` as the Windows entry points.
The Python implementations, shared compiler helper, SDK probe and regression
tests are maintained alongside them. See [BUILD_SDK.md](BUILD_SDK.md) for SDK
setup and recreation on another development machine.

The SDK and wheel scripts share explicit MSVC compiler selection. If an older
wheel run reports GNU 5.3.0 from Anaconda's mingw-w64 directory, rerun the updated
wheel command with the existing SDK. Every wheel attempt uses a fresh temporary
CMake build directory, so no manual cache deletion or SDK rebuild is needed.
Expect the compiler identification to say MSVC. Compiler paths are recorded in
the output directory's build-settings.json.

## One SDK reused across five builds

To build that SDK from scratch, run from a conda-enabled Command Prompt:

```bat
conda activate base
packaging\windows\build_sdk.cmd
```

This creates a root-local **msani_build** conda environment containing Python
3.12, CMake, Ninja, Git, Eigen and native **Boost 1.85.0** from conda-forge.
**RDKit Release_2025_09_1 is compiled from source** with your VS x64 compiler,
without Python wrappers. The old SDK/environment are untouched. An exact conda
package export supports recreation on another development machine; see
[SDK details](BUILD_SDK.md).

The default SDK is installed at:

```text
%USERPROFILE%\msani-sdk\windows-2025.9.1-conda\rdkit
```

After SDK validation passes:

```bat
packaging\windows\build_wheels.cmd --sdk "%USERPROFILE%\msani-sdk\windows-2025.9.1-conda\rdkit" --rdkit-version 2025.9.1 --versions 3.12
```

The generated SDK manifest supplies Boost/Eigen paths automatically. Omit
`--versions 3.12` once it passes to build Python 3.10–3.14.

The wheel script itself does not build an RDKit SDK. Supply your existing Windows shared
SDK, with headers, CMake exports and DLLs, plus matching Boost/Eigen. It accepts
a conda environment root (using its Library directory) or an independent SDK.
Keep this directory between runs: no SDK recompilation is needed per Python.

For an independent SDK, disable RDKit Python wrappers and use consistent MSVC
x64 Release /MD settings. Inspect conda SDK dependencies before assuming they
are Python-independent: version-specific Boost.Python DLLs may prevent reuse.
Linux archives cannot be reused on Windows.

## Run from Windows Terminal's Command Prompt

Activate an environment with Python 3.10 or newer to run the orchestrator:

```bat
conda activate msani
cd /d C:\src\MolSanitizer
packaging\windows\build_wheels.cmd --sdk "C:\msani-sdk\rdkit-2025.9.1" --rdkit-version 2025.9.1 --versions 3.12
```

Replace example paths with your actual checkout and **already built SDK**.
Alternatively, supply the working conda SDK environment:
`--sdk "E:\Anaconda\envs\msani"`. Set `--rdkit-version` to match its actual
C++ SDK version, not merely the Python RDKit package installed alongside it.
This option installs Python RDKit for build checks; it does not update the SDK.

The wrapper initializes VS x64 tools if cl.exe is absent. If that fails, open
**x64 Native Tools Command Prompt for VS**, activate conda there, and retry.
Prefer a local Windows checkout over the WSL network drive for performance.

After 3.12 passes, omit `--versions` to build all five:

```bat
packaging\windows\build_wheels.cmd --sdk "C:\msani-sdk\rdkit-2025.9.1" --rdkit-version 2025.9.1
```

Options:

- `--versions 3.12 3.14`: select a subset for retries.
- `--jobs 2`: compilation concurrency (default 2).
- `--extra-prefixes "C:\boost;C:\eigen"`: additional dependency prefixes.
- `--python-provider launcher`: optional standalone-Python mode, requiring
  installed interpreters and the py launcher, for stricter release tests.

Conda creates fresh temporary environments with only Python and pip, reusing
its normal package-download cache. It never installs another RDKit C++ SDK
into them. Conda uses your configured channels; resolve channel access/terms
or unavailable Python versions before retrying. The SDK and your activated
environment are not modified by the build loop.

## Outputs and tests

Each invocation creates `dist/windows-wheels-<unique>/` containing:

- `wheels/`: repaired wheels, including those whose later tests failed.
- `logs/cp310.log` etc.: build, DLL repair, dependency versions and tests.
- `status.txt`: pass/fail per interpreter.
- `build-settings.json` and `SHA256SUMS`: build inputs and artifact hashes.

Every version gets a fresh CMake directory and build environment. The script
builds a regular wheel, runs delvewheel show/repair and Twine metadata checks,
then installs it into a **separate fresh test environment**. Test PATH contains
only that environment's runtime directories and Windows system directories,
not the original SDK or build environment. Only tests/fixtures are copied;
the source package is not exposed. Native imports are checked before installing
the wheel's pdbqt extra, followed by pip check, dependency listing and full
unittest discovery. Review skips before releasing.

Failures do not prevent trying remaining versions. The script returns nonzero
on failure. Temporary environments are removed; logs and wheels persist.
Nothing is uploaded automatically.

## Publishing gate

Conda is a build tool, not an intended user dependency. The default tests
exclude your build SDK, but conda Python itself supplies runtime libraries.
Passing these tests is therefore not final proof of a conda-free wheel.
Test repaired artifacts with standalone Python on clean Windows 10 before
publishing, with no developer SDK directories available.

Review DLL license notices and version-specific Python/Boost.Python
dependencies. The win_amd64 tag does not encode a minimum Windows version.
ARM64, 32-bit and free-threaded Python are outside this matrix.

## Development checks

The packaging helper tests do not need Windows, a compiler, or network access:

```sh
python -m unittest discover -s packaging/windows -p 'test_*.py' -v
```

Keep these scripts, tests and the C++ SDK probe in Git. Generated wheels, logs,
build directories and Python bytecode are already excluded by the repository's
ignore rules. SDKs live outside the checkout and should not be committed.

These builders target Windows x64 only. Apple Silicon requires a separate
macOS arm64 SDK and wheel build; Windows DLLs and win_amd64 wheels cannot be
used on macOS. A macOS packaging workflow is not provided by this directory.
