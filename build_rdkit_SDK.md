# MolSanitizer — Publishing to PyPI: Complete Guide

> [!NOTE]
> This guide assumes you're working in WSL with the `msani` conda environment.
> All commands use `source /home/phonglam/anaconda3/etc/profile.d/conda.sh && conda activate msani` as the setup step.

---

## Table of Contents

1. [Build Static RDKit SDK (one-time)](#1-build-static-rdkit-sdk-one-time)
2. [Build the Self-Contained Wheel](#2-build-the-self-contained-wheel)
3. [Verify the Wheel](#3-verify-the-wheel)
4. [Set Up PyPI Accounts](#4-set-up-pypi-accounts)
5. [Upload to TestPyPI](#5-upload-to-testpypi)
6. [Upload to PyPI (Production)](#6-upload-to-pypi-production)
7. [Automate with GitHub Actions](#7-automate-with-github-actions)
8. [Checklist for Each Release](#8-checklist-for-each-release)

---

## 1. Build Static RDKit SDK (One-Time)

This builds RDKit 2025.09.5 as static `.a` libraries so they get baked into your wheel.
You only need to do this once (or when upgrading RDKit version).

> [!NOTE]
> On this workstation the SDK is already installed at
> `/home/phonglam/msani-rdkit-sdk`. It has been checked for RDKit 2025.09.5
> headers, CMake metadata, and static archives. Do not rebuild it unless the
> RDKit version or compiler environment changes.

```bash
# ── Prerequisites ──
# Make sure you have build tools (these are already in msani env)
# cmake, ninja-build, eigen3, boost headers, g++

# ── Choose where to install the SDK (your home dir, no sudo needed) ──
export RDKIT_SDK=$HOME/msani-rdkit-sdk

# ── Clone the pinned RDKit version into /tmp (temporary, just for building) ──
cd /tmp
git clone --depth 1 --branch Release_2025_09_5 \
    https://github.com/rdkit/rdkit.git rdkit-src

# ── Configure: build in /tmp, but CMAKE_INSTALL_PREFIX tells cmake
#    where to COPY the final headers + libraries when you run "cmake --install" ──
cmake -S rdkit-src -B rdkit-build -G Ninja \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX=$RDKIT_SDK \
    -DRDK_INSTALL_INTREE=OFF \
    -DRDK_BUILD_STATIC_LIBS_ONLY=ON \
    -DRDK_BUILD_PYTHON_WRAPPERS=OFF \
    -DRDK_BUILD_CPP_TESTS=OFF \
    -DRDK_BUILD_INCHI_SUPPORT=OFF \
    -DCMAKE_POSITION_INDEPENDENT_CODE=ON

# ── Compile (uses all cores; takes ~10-20 minutes) ──
cmake --build rdkit-build --parallel $(nproc)

# ── Install: copies the built .a files and headers from /tmp/rdkit-build
#    into ~/msani-rdkit-sdk/  (no sudo needed since it's in your home dir) ──
cmake --install rdkit-build

# ── Clean up the build tree (optional, frees ~2 GB) ──
rm -rf /tmp/rdkit-src /tmp/rdkit-build
```

> [!IMPORTANT]
> `-DCMAKE_POSITION_INDEPENDENT_CODE=ON` is critical! Static archives will be linked into
> shared `.so` extension modules, which requires PIC (Position Independent Code).

After this you'll have:
```
~/msani-rdkit-sdk/
├── include/rdkit/       ← headers
├── lib/
│   ├── libRDKitGraphMol.a
│   ├── libRDKitForceField.a
│   └── ...              ← all static archives
└── lib/cmake/rdkit/     ← CMake config files
```

RDKit 2025.09.5 actually names these archives
`libRDKitGraphMol_static.a`, `libRDKitForceField_static.a`, and so on, and
exports matching CMake targets such as `RDKit::GraphMol_static`.

---

## 2. Build the Self-Contained Wheel

```bash
cd /home/phonglam/github/MolSanitizer

# Build the wheel with static linkage. This form works with the tools already
# installed in the msani environment and does not download build dependencies.
mkdir -p dist/local
RDBASE=/home/phonglam/msani-rdkit-sdk \
CMAKE_ARGS="-DMSANI_RDKIT_LINKAGE=static -DCMAKE_PREFIX_PATH=/home/phonglam/msani-rdkit-sdk\;/home/phonglam/anaconda3/envs/msani" \
python -m pip wheel . --no-build-isolation --no-deps --wheel-dir dist/local

# The wheel lands in dist/
ls dist/local/
# → molsanitizer-0.7.0-cp312-cp312-linux_x86_64.whl
```

### Verify No Dynamic RDKit Dependencies

```bash
# Check what the .so files inside the wheel link to
unzip -o dist/local/molsanitizer-0.7.0-cp312-cp312-linux_x86_64.whl \
    -d /tmp/wheel-check

# Should show NO libRDKit*.so dependencies:
ldd /tmp/wheel-check/msani_confgen_cpp*.so 2>&1 | grep -i rdkit
ldd /tmp/wheel-check/msani_stereoisomers*.so 2>&1 | grep -i rdkit
ldd /tmp/wheel-check/amsolcpp/_amsolcpp*.so 2>&1 | grep -i rdkit

# If these commands print nothing → success! RDKit is statically linked.
# If they still show libRDKit*.so → something went wrong, re-check the build.
```

### Run `auditwheel` (for PyPI-uploadable wheels)

`auditwheel repair` also needs the `patchelf` executable. Install it in the
build environment first if `command -v patchelf` prints nothing:

```bash
conda install -c conda-forge patchelf
```

```bash
# Check the wheel's compatibility
auditwheel show dist/local/molsanitizer-0.7.0-cp312-cp312-linux_x86_64.whl

# Repair it to make it manylinux-compliant
# (vendors any remaining non-standard shared libs into the wheel)
mkdir -p dist/repaired
auditwheel repair dist/local/molsanitizer-0.7.0-cp312-cp312-linux_x86_64.whl \
    --plat manylinux_2_39_x86_64 \
    -w dist/repaired/

# The repaired wheel is in dist/repaired/
ls dist/repaired/
```

> [!TIP]
> If `auditwheel` says the wheel is already compliant, the original wheel is fine to upload.
> If it finds vendored libs that need bundling, the repaired wheel in `dist/repaired/` is what you upload.

> [!CAUTION]
> `manylinux_2_39` is suitable for the current TestPyPI prototype, but it will
> not install on older Linux distributions. The production wheel must be built
> from `packaging/linux/Dockerfile` in a manylinux 2.28 container; changing the
> tag alone cannot make a WSL-built binary compatible with older glibc.

### Supporting Ubuntu 20.04

Ubuntu 20.04 uses glibc 2.31. A wheel built and repaired as
`manylinux_2_28_x86_64` requires glibc 2.28 or newer, so it supports Ubuntu
20.04 as well as newer compatible Linux distributions. The current local
`manylinux_2_39_x86_64` prototype does **not** support Ubuntu 20.04.

Build the release wheel inside `quay.io/pypa/manylinux_2_28_x86_64`, as shown
in the GitHub Actions section. Do not relabel a locally built 2.39 wheel as
2.28; all native code, including the static RDKit SDK, must be compiled inside
the manylinux 2.28 environment.

The Python tag is independent of Linux compatibility. A `cp312-cp312` wheel
requires CPython 3.12, including on Ubuntu 20.04. Build separate `cp310`,
`cp311`, and `cp312` wheels if you want to support all Python versions allowed
by `pyproject.toml`. Ubuntu's system Python does not determine which wheel a
user may install; the interpreter running `pip` does.

---

## 3. Verify the Wheel

**This is the most important step.** Test in a completely clean environment with no conda:

```bash
# Leave both an old test venv and Conda before selecting the interpreter.
deactivate 2>/dev/null || true
conda deactivate 2>/dev/null || true

# Ubuntu/WSL needs this package once. If the venv command below reports that
# ensurepip is unavailable, run:
# sudo apt update && sudo apt install python3.12-venv

# Recreate the test environment explicitly from Ubuntu's Python, not whichever
# python3.12 happens to be first on PATH in Conda.
rm -rf /tmp/test-msani-wheel
/usr/bin/python3.12 -m venv /tmp/test-msani-wheel
source /tmp/test-msani-wheel/bin/activate

# This must print /usr. A Conda path means the venv is not a clean test.
python -c "import sys; print('base Python:', sys.base_prefix); assert sys.base_prefix == '/usr'"

# Install ONLY from the wheel + PyPI dependencies
python -m pip install \
    /home/phonglam/github/MolSanitizer/dist/repaired/molsanitizer-*.whl

# Run outside the repository so its source tree cannot shadow the installed
# package. This calls the native MolSanitizer conformer implementation.
cd /tmp
python -c "
import msani_confgen_cpp
import msani_stereoisomers
from amsolcpp import _amsolcpp
from rdkit import Chem, rdBase
print('rdkit version:', rdBase.rdkitVersion)

# Test conformer generation (exercises confgen C++ module)
mol = Chem.MolFromSmiles('CCO')
mol = Chem.AddHs(mol)
result = msani_confgen_cpp.embed_multiple_confs(mol, 2, 'etkdgv3')
assert result.GetNumConformers() > 0
print('Native conformers:', result.GetNumConformers())
print('MolSanitizer loaded successfully')
"

# Run the full test suite while still outside the repository. Test data are
# resolved relative to test_msani.py, but imports come from the installed wheel.
python -m unittest discover \
    -s /home/phonglam/github/MolSanitizer/test \
    -p 'test_msani.py' \
    -v

deactivate
rm -rf /tmp/test-msani-wheel
```

> [!CAUTION]
> If the clean venv test fails with `ImportError: libRDKit*.so not found`, the static linkage
> didn't work. Go back to step 2 and make sure `MSANI_RDKIT_LINKAGE=static` was set and
> `CMAKE_PREFIX_PATH` pointed to your static SDK.

> [!CAUTION]
> If `pyvenv.cfg` names `/home/phonglam/anaconda3/...` as `home` or
> `sys.base_prefix`, delete and recreate the venv with `/usr/bin/python3.12`.
> Mixing the Ubuntu executable with Conda's `_ctypes` module produces errors
> such as `undefined symbol: _PyErr_SetLocaleString`.

---

## 4. Set Up PyPI Accounts

### Create Accounts

1. **TestPyPI** (sandbox): Go to [https://test.pypi.org/account/register/](https://test.pypi.org/account/register/)
2. **PyPI** (production): Go to [https://pypi.org/account/register/](https://pypi.org/account/register/)

### Create API Tokens

For **both** sites:

1. Log in → Go to **Account settings** → **API tokens**
2. Click **"Add API token"**
3. Token name: `MolSanitizer-upload` (or anything descriptive)
4. Scope: **Entire account** (for first upload; after that you can scope to the project)
5. **Copy the token immediately** — you won't see it again!

Tokens look like: `pypi-AgEIcHlwaS5vcmcC...` (very long string)

### Store Tokens Securely

Create a `~/.pypirc` file:

```ini
[distutils]
index-servers =
    pypi
    testpypi

[pypi]
username = __token__
password = pypi-AgEIcHlwaS5vcmcC...YOUR_REAL_PYPI_TOKEN...

[testpypi]
repository = https://test.pypi.org/legacy/
username = __token__
password = pypi-AgEIcHlwaS5vcmcC...YOUR_REAL_TESTPYPI_TOKEN...
```

```bash
chmod 600 ~/.pypirc   # Restrict permissions
```

### Install `twine`

```bash
pip install twine
```

---

## 5. Upload to TestPyPI

**Always test here first!** TestPyPI is a sandbox — nothing you upload here affects real users.

```bash
# Check the wheel is well-formed
twine check dist/repaired/MolSanitizer-*.whl

# Upload to TestPyPI
twine upload --repository testpypi dist/repaired/MolSanitizer-*.whl
```

### Test the TestPyPI Upload

```bash
# Create another clean venv
python3.12 -m venv /tmp/test-testpypi
source /tmp/test-testpypi/bin/activate

# Install from TestPyPI, but get dependencies from real PyPI
pip install --index-url https://test.pypi.org/simple/ \
    --extra-index-url https://pypi.org/simple/ \
    MolSanitizer==0.7.0

python -c "import msani; print(msani.__version__)"

deactivate
rm -rf /tmp/test-testpypi
```

> [!NOTE]
> TestPyPI often doesn't have all dependencies (numpy, pandas, etc.), so `--extra-index-url`
> tells pip to fall back to real PyPI for those. This is normal and expected.

### Common Issues

| Problem | Solution |
|---------|----------|
| `403 Forbidden` | Token is wrong or expired — regenerate it |
| `400 File already exists` | You can't re-upload the same version. Bump the version number (e.g., `0.7.0` → `0.7.1`) |
| `Invalid distribution` | Run `twine check` and fix any warnings |

---

## 6. Upload to PyPI (Production)

Once TestPyPI works perfectly:

```bash
# Final check
twine check dist/repaired/MolSanitizer-*.whl

# Upload to REAL PyPI
twine upload dist/repaired/MolSanitizer-*.whl
```

Your package is now live at `https://pypi.org/project/MolSanitizer/`!

Users can install it with:
```bash
pip install MolSanitizer
```

> [!CAUTION]
> **You cannot delete or overwrite a version on PyPI once uploaded.**
> If you find a bug, you must release a new version (e.g., `0.7.1`).
> This is why TestPyPI exists — always test there first!

---

## 7. Automate with GitHub Actions

Once manual publishing works, automate it. Create `.github/workflows/publish.yml`:

```yaml
name: Build & Publish to PyPI

on:
  release:
    types: [published]     # Triggers when you create a GitHub Release

  workflow_dispatch:        # Manual trigger for testing
    inputs:
      upload_target:
        description: "Upload target"
        required: true
        default: "testpypi"
        type: choice
        options:
          - testpypi
          - pypi
          - none

jobs:
  build-wheel:
    runs-on: ubuntu-latest
    container:
      image: quay.io/pypa/manylinux_2_28_x86_64

    steps:
      - uses: actions/checkout@v4

      # ── Build static RDKit SDK ──
      - name: Install build dependencies
        run: dnf -y install git cmake ninja-build boost-devel eigen3-devel

      - name: Build static RDKit SDK
        run: |
          git clone --depth 1 --branch Release_2025_09_5 \
              https://github.com/rdkit/rdkit.git /tmp/rdkit
          cmake -S /tmp/rdkit -B /tmp/rdkit-build -G Ninja \
              -DCMAKE_BUILD_TYPE=Release \
              -DCMAKE_INSTALL_PREFIX=/opt/msani-rdkit \
              -DRDK_INSTALL_INTREE=OFF \
              -DRDK_BUILD_STATIC_LIBS_ONLY=ON \
              -DRDK_BUILD_PYTHON_WRAPPERS=OFF \
              -DRDK_BUILD_CPP_TESTS=OFF \
              -DRDK_BUILD_INCHI_SUPPORT=OFF \
              -DCMAKE_POSITION_INDEPENDENT_CODE=ON
          cmake --build /tmp/rdkit-build --parallel $(nproc)
          cmake --install /tmp/rdkit-build

      # ── Build the wheel ──
      - name: Build wheel
        run: |
          /opt/python/cp312-cp312/bin/pip install build scikit-build-core pybind11
          /opt/python/cp312-cp312/bin/python -m build --wheel \
              -C cmake.define.CMAKE_PREFIX_PATH=/opt/msani-rdkit \
              -C cmake.define.MSANI_RDKIT_LINKAGE=static

      - name: Repair wheel with auditwheel
        run: |
          /opt/python/cp312-cp312/bin/pip install auditwheel
          /opt/python/cp312-cp312/bin/auditwheel repair dist/*.whl \
              --plat manylinux_2_28_x86_64 -w dist/repaired/

      - name: Upload wheel artifact
        uses: actions/upload-artifact@v4
        with:
          name: wheel-linux-cp312
          path: dist/repaired/*.whl

  # ── Publish ──
  publish:
    needs: build-wheel
    runs-on: ubuntu-latest
    permissions:
      id-token: write   # Required for trusted publishing (no token needed!)

    steps:
      - uses: actions/download-artifact@v4
        with:
          name: wheel-linux-cp312
          path: dist/

      # For GitHub Release → publish to PyPI
      - name: Publish to PyPI
        if: github.event_name == 'release'
        uses: pypa/gh-action-pypi-publish@release/v1
        with:
          packages-dir: dist/

      # For manual trigger → publish to chosen target
      - name: Publish to TestPyPI
        if: >-
          github.event_name == 'workflow_dispatch' &&
          github.event.inputs.upload_target == 'testpypi'
        uses: pypa/gh-action-pypi-publish@release/v1
        with:
          repository-url: https://test.pypi.org/legacy/
          packages-dir: dist/
```

### Set Up Trusted Publishing (Recommended over tokens)

Instead of managing API tokens, PyPI supports **Trusted Publishing** via GitHub Actions:

1. Go to [pypi.org](https://pypi.org) → Your project → **Publishing** → **Add a new publisher**
2. Fill in:
   - **Owner**: `phonglam3103`
   - **Repository**: `MolSanitizer`
   - **Workflow name**: `publish.yml`
   - **Environment**: _(leave blank)_
3. Do the same on [test.pypi.org](https://test.pypi.org) for TestPyPI

Now the GitHub Action can publish **without any tokens** — GitHub's OIDC proves it's your repo.

---

## 8. Checklist for Each Release

```
Pre-release:
  [ ] Update version in pyproject.toml (e.g., 0.7.0 → 0.8.0)
  [ ] Update CHANGELOG.md
  [ ] All tests pass: python -m pytest test/ -v
  [ ] Commit everything, push to main

Build:
  [ ] Build static wheel (step 2)
  [ ] auditwheel repair (step 2)
  [ ] Verify in clean venv (step 3)

Publish:
  [ ] Upload to TestPyPI first (step 5)
  [ ] Test install from TestPyPI
  [ ] Upload to PyPI (step 6) — or create GitHub Release to trigger CI

Post-release:
  [ ] Verify: pip install MolSanitizer==X.Y.Z in a clean env
  [ ] Tag the release: git tag v0.8.0 && git push --tags
  [ ] Create GitHub Release from the tag
```

---

## Quick Reference: Key Commands

```bash
# ── Build ──
RDBASE=/home/phonglam/msani-rdkit-sdk \
CMAKE_ARGS="-DMSANI_RDKIT_LINKAGE=static -DCMAKE_PREFIX_PATH=/home/phonglam/msani-rdkit-sdk\;/home/phonglam/anaconda3/envs/msani" \
python -m pip wheel . --no-build-isolation --no-deps --wheel-dir dist/local

# ── Check ──
auditwheel show dist/local/*.whl
twine check dist/repaired/*.whl

# ── Repair ──
auditwheel repair dist/local/*.whl --plat manylinux_2_39_x86_64 -w dist/repaired/

# ── Upload ──
twine upload --repository testpypi dist/repaired/*.whl    # TestPyPI
twine upload dist/repaired/*.whl                           # Real PyPI
```
