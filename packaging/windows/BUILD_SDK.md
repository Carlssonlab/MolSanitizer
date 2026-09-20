# Reusable Windows RDKit SDK

Build RDKit **2025.9.1** from source with MSVC, using native Boost **1.85.0**
libraries and headers from conda-forge. This follows the dependency strategy in
[RDKit's Windows CI](https://github.com/rdkit/rdkit/blob/master/.azure-pipelines/vs_build.yml).
There is no Boost source build or b2 setup step.

## First build

On a Windows x64 development machine, install conda and Visual Studio 2022 Build
Tools with the C++ workload and Windows SDK. Open a conda-enabled **Command
Prompt** in the MolSanitizer checkout. The launching Python must be 3.10 or newer.

```bat
conda activate base
packaging\windows\build_sdk.cmd
```

The wrapper initializes x64 MSVC. If that fails, use an x64 Native Tools Command
Prompt and activate conda there. Docker is not required. Internet access and
several GB of writable local disk are needed.

New default root: `%USERPROFILE%\msani-sdk\windows-2025.9.1-conda`.

The old SDK root and globally named msani_build environment are left untouched.
The new environment lives at `<new-root>\msani_build`. You need not activate it
manually. Do not modify it after building.

```bat
packaging\windows\build_sdk.cmd --root C:\msani-sdk\windows-2025.9.1-conda --jobs 4
```

## What the script does

1. Creates/reuses the root-local conda environment with libboost=1.85.0,
   libboost-devel=1.85.0, Python 3.12, Eigen 3.4, CMake 3.28–3.x, Ninja and Git.
   It does not install Python RDKit or Boost.Python into the SDK environment.
2. Exports exact resolved conda package URLs/checksums to conda-win-64.lock.txt,
   and records compiler and Windows SDK versions in toolchain.json.
3. Clones Release_2025_09_1 and records its commit. Builds shared RDKit with
   MSVC, Ninja, Release and /MD, Python/SWIG wrappers disabled, and thread-safe
   support enabled. The SDK has a separate install prefix.
4. Runs a C++ consumer probe and checks RDKit/Boost DLL imports for Python
   dependencies. This is a lightweight SDK check, not RDKit's full CTest suite.
5. Copies RDKit and cached Boost license texts, then writes
   rdkit/msani-sdk.json only after validation succeeds.

Conda Boost can introduce compression/runtime dependencies even though RDKit's
compressed suppliers are disabled. Wheel repair must resolve transitive DLLs;
do not copy just the RDKit DLLs manually. Do not run conda clean before SDK
completion: Boost notices are collected from the extracted package cache.

## Build and test wheels

After the SDK reports success:

```bat
packaging\windows\build_wheels.cmd --sdk "%USERPROFILE%\msani-sdk\windows-2025.9.1-conda\rdkit" --rdkit-version 2025.9.1 --versions 3.12
```

Omit `--versions 3.12` to build Python 3.10–3.14. The same SDK serves every
Python version. The wheel script reads the manifest to locate native dependencies,
repairs the wheels with delvewheel, and runs unit tests.
Only distribute repaired, tested wheels. End users install those using pip;
they do not need this SDK, conda, or a compiler.

## Recreate on another development machine

Save these with your release build records:

- The MolSanitizer repository commit used for the build.
- conda-win-64.lock.txt (exact conda package artifacts).
- toolchain.json (MSVC toolset and Windows SDK versions).
- rdkit/msani-sdk.json (including the RDKit commit).

On the second Windows x64 machine, install conda and matching Visual Studio
toolset/Windows SDK versions. Check out the same MolSanitizer commit. Copy the
lock file outside the new SDK root, then run:

```bat
conda activate base
packaging\windows\build_sdk.cmd --root C:\msani-sdk\windows-2025.9.1-conda --conda-lock C:\build-records\conda-win-64.lock.txt --rdkit-commit RECORDED_FULL_COMMIT_SHA
```

Replace RECORDED_FULL_COMMIT_SHA with rdkit_commit from the saved manifest.
Conda recreates the exact package artifacts; the script rejects a different
RDKit commit. Paths may differ between machines. Rebuild locally instead of
copying CMake caches or an installed SDK containing absolute paths.

Without --conda-lock, conda may choose newer builds of dependencies.
The lock does not install Visual Studio or guarantee byte-identical binaries.
RDKit can download auxiliary sources during configuration: this is not an offline
or fully hermetic build. Preserve successful sources/package caches if long-term
availability matters.

## Cache and release checks

The bootstrap explicitly pins both C and C++ to MSVC cl.exe; Ninja alone does
not select MSVC. Inherited GCC compiler variables/flags are cleared. Builds use
rdkit-build-msvc and probe-build-msvc so an earlier MinGW CMake cache is not
reused. If an earlier run selected Anaconda's mingw-w64 compiler, simply rerun
the updated command with the same SDK root. The environment and downloaded
RDKit source are retained, and the old build directories are left untouched.

Rerunning with the same root reuses the environment, RDKit source and incremental
build. Compiler/environment changes are rejected; choose a new root instead.
No existing SDK directories are automatically deleted. Concurrent builds sharing
one root are not supported.

SDK validation alone does not establish Windows 10 compatibility. Before release,
test repaired wheels on clean target Windows machines without conda or the SDK,
and review notices for all bundled DLLs, including transitive dependencies.
