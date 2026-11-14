# Releasing MolSanitizer wheels

**Prepare release** opens one changelog PR. Once it is merged, pushing a
bare release tag starts **Wheels**, which builds, verifies, and publishes the
release to PyPI using Trusted Publishing. The **SDK** workflow provisions
build inputs separately.

## Development and release versions

Git tags are the source of the package version. `setuptools-scm` resolves the
version through the `scikit-build-core` metadata provider; there is no static
version to bump in `pyproject.toml` or `docs/source/conf.py`.

| Checkout | Example version |
| --- | --- |
| Clean commit tagged `0.7.5` | `0.7.5` |
| Later commits on `main` | `0.7.6.devN+g<commit>` |
| Uncommitted changes after those commits | `0.7.6.devN+g<commit>.dYYYYMMDD` |

`N` is the distance from the release tag. A dirty checkout at the tag itself
is also a development build, such as `0.7.6.dev0+g<commit>.dYYYYMMDD`. SDK tags
are excluded from version detection. Release tags use bare version numbers, without a `v`
prefix; release candidates such as `0.7.5rc1` are also supported. Development
versions infer the next patch version; the final tag can choose a different
minor or major version when appropriate.

To release `0.7.5`:

1. Run **Prepare release** with `0.7.5` and review its changelog PR.
2. Merge that PR, then tag its exact merge commit:

   ```bash
   git tag 0.7.5 <release-merge-commit-sha>
   git push origin 0.7.5
   ```

3. Wait for **Wheels** to pass and publish. Read the Docs can build `stable`
   from that tag. Later commits automatically get development versions, so
   there is no second version-bump PR. While `main` still points to the tagged
   commit, it correctly reports the release version.

Use a full clone with release tags for local source builds. With the build
requirements installed, `python -m setuptools_scm` prints the resolved version.
Source distributions preserve it in `PKG-INFO`; Git archives preserve it via
`.git_archival.txt` and `.gitattributes`. An arbitrary source copy without Git,
archive metadata, or `PKG-INFO` fails version detection instead of reporting an
unrelated release. The Linux wheel helper resolves the original version before
copying source without `.git`. Existing installations need reinstalling to
pick up changed package metadata.

## Read the Docs versions

Use `latest` for development documentation built from `main`, and `stable`
for documentation built from the most recent final release tag. Sphinx uses
the same Git version configuration as the package, so `latest` can show
`0.7.6.devN+g<commit>` while `stable` shows `0.7.5`. The Read the Docs build
hook fetches full history and tags before resolving this version; no native
extension build is required to determine the documentation version.

Configure this in the Read the Docs project dashboard:

1. Check that the repository URL points to the repository where release tags
   are pushed, and set **Default branch** to `main` in the project settings.
2. In **Versions**, activate the latest final release tag, `stable`, and
   `latest`, and leave them visible (not hidden). Check the Git reference
   behind `stable` and ensure its build succeeds. Trigger a build if newly
   pushed tags have not been synced. A missing release tag must point to the
   actual published release commit, not the current development commit.
3. Set **Default version** to `stable` in the project settings. This makes
   the documentation root URL open the released documentation.
4. Under **Settings > Addons > Notifications**, enable the **Latest version
   warning** to identify development documentation to readers.
5. Under **Automation Rules**, add a rule with **Custom match**
   `^[0-9]+\.[0-9]+\.[0-9]+$`, **Version type** `Tag`, and **Action**
   `Activate version`. This activates future final release tags and excludes
   SDK tags, release candidates, and development versions. Existing versions
   must be activated manually.
6. Verify that `/en/stable/` shows the released version, `/en/latest/` shows
   the development version after commits beyond the release tag are merged,
   and the root URL redirects to `/en/stable/`. Check that the README's stable guide links
   resolve; older release tags may have different documentation pages.

These project settings are managed in the dashboard, separately from the
build configuration in `.readthedocs.yaml`. Future commits on `main`
do not change the documentation stored at existing release tags. Read the
Docs selects `stable` from final version tags, independently of whether a
package upload succeeded; coordinate release tags with publication.

The sidebar shows the Sphinx release number and adds it to the selected
`stable` or `latest` option, for example `0.7.6.dev12 (latest)`. The label below
and the selector tooltip retain the full version including its commit suffix.
The release number remains visible when the selector or JavaScript is unavailable.
After merging documentation styling changes, rebuild `latest`. The `stable`
build receives those changes when its selected release contains them; rebuilding
an older release tag continues to use that tag's templates and styles.

See the Read the Docs guides for
[versions](https://docs.readthedocs.com/platform/latest/versions.html),
[automation rules](https://docs.readthedocs.com/platform/stable/guides/automation-rules.html),
and [notifications](https://docs.readthedocs.com/platform/latest/doc-notifications.html).

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

The **RDKit compatibility** workflow (`rdkit-compatibility.yml`) runs on code
PRs, pushes to `main`, and manual dispatch. It restores the existing
`sdk-2026.09.2` Linux SDK and builds one repaired Linux x86_64 CPython 3.12
wheel against native RDKit `2026.9.2`. At the start of each run, a discovery
job queries PyPI's JSON simple index and selects every stable RDKit release
from `2025.3.1` onward with a non-yanked wheel compatible with Linux x86_64
CPython 3.12. It checks wheel tags and `Requires-Python` inside the same
`python:3.12-slim` image used for testing. Alpha, beta, release-candidate,
development, local, and source-only releases are excluded; stable postreleases
are included. New compatible releases are added automatically on the next run.

The discovery job records the exact version list in its job summary and
passes it to the test matrix. One job per version downloads the same
MolSanitizer artifact and runs the full installed-wheel suite. Each runtime
is constrained during both core and PDBQT extras
installation, then checked against the imported RDKit module and package
metadata before running the suite. Tests run from a temporary directory using
installed extensions. Failures in one runtime do not cancel the other cases.
RDKit is installed from binary wheels only. No SDK or MolSanitizer rebuild
occurs between runtime versions. A PyPI lookup failure or empty selection
fails the workflow instead of silently dropping version coverage. The total
runner time grows with the number of compatible RDKit releases; parallelism
depends on the repository's available GitHub Actions runners.

The workflow requires the SDK release and its `msani-sdk-linux-x86_64.tar.gz`
asset in the repository running the workflow. It does not publish packages.
The existing **CI** workflow continues its baseline platform checks. To check
the runner's version constraints and failure handling locally without native
builds or downloads, run:

```bash
python -m pip install 'packaging>=23'
python -m unittest discover -s packaging/linux -p 'test_*.py'
```

To print the live version list locally, run this in a Linux x86_64 CPython
3.12 environment (ideally the same `python:3.12-slim` image):

```bash
python packaging/linux/discover_rdkit_versions.py --minimum 2025.3.1
```

The **Wheels** workflow builds, repairs and clean-install tests every wheel.

| Trigger | What runs |
| --- | --- |
| `workflow_dispatch`, scope `linux-only` | Linux x86_64, CPython 3.12 |
| `workflow_dispatch`, scope `full` | 25 jobs: 5 targets x CPython 3.10-3.14 |
| Version tag, e.g. `0.7.5` | Full matrix, plus a tag/version consistency check |

Keeping PR compatibility checks on Linux is deliberate. A full matrix costs roughly $7-10 of the
monthly allowance, most of it macOS; a per-commit full matrix would exhaust
3,000 included minutes in under two runs.

Every wheel is repaired (auditwheel / delocate / delvewheel) and then
installed into a clean interpreter that has no SDK, conda or build directory
on its search path, where the repository suite runs against it. On Linux the
test runs in a bare `python:X.Y-slim` image rather than the manylinux build
image, because a wheel needing build-image libraries is not redistributable.

## 3. Publish the verified bundle

A release tag push publishes the wheels in the `molsanitizer-dist` bundle to PyPI only after
all builds, clean-install tests, and version checks pass. Configure PyPI
Trusted Publishing for `.github/workflows/wheels.yml` and the protected
GitHub environment `pypi`. Review the workflow's publishing job and its result
before treating the release as published. The sdist remains in the download
bundle; it is excluded from automatic uploads because it needs an RDKit C++ SDK.

For TestPyPI, run **Wheels** manually from a release tag with scope `full` and
`publish_testpypi` enabled. Configure its Trusted Publisher and the `testpypi`
GitHub environment. Manual runs do not trigger the PyPI publishing job.
Select a release tag rather than `main`: development versions include a local
commit suffix, which package indexes do not accept.

The downloadable bundle can also be checked locally:

```bash
unzip molsanitizer-dist.zip -d dist
python -m pip install --upgrade twine
python -m twine check --strict dist/*
(cd dist && sha256sum -c SHA256SUMS)
python packaging/check_release_versions.py dist --expected-version 0.7.5
```

Install from TestPyPI on a clean machine to check the uploaded release:

```bash
python -m pip install --index-url https://test.pypi.org/simple/ \
  --extra-index-url https://pypi.org/simple/ MolSanitizer==0.7.5
```

## Notes

- **Claim the name first.** PyPI has no reservation mechanism: a project name
  is claimed by uploading a distribution. Upload to TestPyPI and PyPI early,
  while the repository is still private.
- **The sdist is uploadable but not installable** without an RDKit C++ SDK.
  Publishing it means `pip install` on an unsupported platform attempts a
  source build and fails. Omit `dist/*.tar.gz` from the upload if that
  trade-off is not wanted.
- **Versions come from tags.** The collect job checks each wheel's `METADATA`
  and the source distribution's `PKG-INFO`, requiring every artifact to match
  the release tag. Manual bundles must also agree on a single version.
- **Full Git history is required.** Build jobs and Read the Docs fetch release
  tags so development versions are computed consistently.
