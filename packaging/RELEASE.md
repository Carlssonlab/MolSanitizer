# Releasing MolSanitizer wheels

Two workflows, and one manual upload. Nothing is published automatically:
GitHub Actions only produces a verified bundle, and a person uploads it.

## 1. Build the RDKit SDK (rarely)

The RDKit C++ SDK is an input to every wheel, never part of one. It is built
once per pin and attached to a permanent pre-release, because at macOS runner
rates an SDK that has to be recompiled per run costs more per month than the
entire wheel matrix.

Run the **SDK** workflow (`workflow_dispatch`) when, and only when:

- the RDKit pin changes (`RDKIT_REF` / `rdkit_version`),
- Boost, Eigen or another SDK dependency changes,
- `packaging/linux/build_sdk.sh`, `packaging/macos/build_sdk.sh` or
  `packaging/windows/build_sdk.py` changes.

Pick `all` for a new pin, or a single target to repair one asset. The run
creates the `sdk-<version>` pre-release if it does not exist and attaches:

```
msani-sdk-linux-x86_64.tar.gz     static RDKit + Boost + Eigen under /opt
msani-sdk-linux-aarch64.tar.gz    static RDKit + Boost + Eigen under /opt
msani-sdk-macos-arm64.tar.gz      conda-forge prefix, workspace-relative
msani-sdk-macos-x86_64.tar.gz     conda-forge prefix, workspace-relative
msani-sdk-windows-x64.tar.gz      shared MSVC SDK + Boost, C:\msani-sdk\...
```

Two of these are position-dependent by nature. A conda prefix is not
relocatable, so the macOS assets unpack at the same workspace-relative path
they were built at; `packaging/windows/build_sdk.py` records absolute
dependency prefixes in `msani-sdk.json`, so the Windows asset unpacks at a
fixed `C:\msani-sdk\` path. Both hold on hosted runners, where those paths are
stable. **These assets are CI inputs; for local builds, run the platform's own
`build_sdk.sh` / `build_sdk.cmd`, which provision the SDK in place.**

## 2. Build the wheels

The **Wheels** workflow builds, repairs and clean-install tests every wheel.

| Trigger | What runs |
| --- | --- |
| Pull request | One job: Linux x86_64, CPython 3.12 (trigger enabled once an SDK release exists) |
| `workflow_dispatch`, scope `linux-only` | Same cheap path |
| `workflow_dispatch`, scope `full` | 25 jobs: 5 targets x CPython 3.10-3.14 |
| Tag `v*` | Full matrix, plus a tag/version consistency check |

The cheap default is deliberate. A full matrix costs roughly $7-10 of the
monthly allowance, most of it macOS; a per-commit full matrix would exhaust
3,000 included minutes in under two runs.

Every wheel is repaired (auditwheel / delocate / delvewheel) and then
installed into a clean interpreter that has no SDK, conda or build directory
on its search path, where the repository suite runs against it. On Linux the
test runs in a bare `python:X.Y-slim` image rather than the manylinux build
image, because a wheel needing build-image libraries is not redistributable.

## 3. Upload by hand

Download the `molsanitizer-dist` artifact from the run summary, then:

```bash
unzip molsanitizer-dist.zip -d dist
python -m pip install --upgrade twine
python -m twine check --strict dist/*
sha256sum -c dist/SHA256SUMS

python -m twine upload --repository testpypi dist/*
```

Install from TestPyPI on a clean machine before touching real PyPI:

```bash
python -m pip install --index-url https://test.pypi.org/simple/ \
  --extra-index-url https://pypi.org/simple/ MolSanitizer
```

Only then:

```bash
python -m twine upload dist/*
```

## Notes

- **Claim the name first.** PyPI has no reservation mechanism: a project name
  is claimed by uploading a distribution. Upload to TestPyPI and PyPI early,
  while the repository is still private.
- **The sdist is uploadable but not installable** without an RDKit C++ SDK.
  Publishing it means `pip install` on an unsupported platform attempts a
  source build and fails. Omit `dist/*.tar.gz` from the upload if that
  trade-off is not wanted.
- **Version bumps** are a `pyproject.toml` edit. A `v*` tag whose name does
  not match `project.version` fails the collect job before any bundle is
  treated as releasable.
- **Switching to Trusted Publishing** once the repository is public needs a
  publish job with `permissions: id-token: write`, a protected `pypi`
  environment and `pypa/gh-action-pypi-publish`. Nothing else in these
  workflows has to change.
