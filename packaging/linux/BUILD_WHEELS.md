# Build and test Linux wheels locally

Run from Ubuntu/WSL on your Docker host, **not inside a container**:

```bash
cd /home/phonglam/github/MolSanitizer
bash packaging/linux/build_wheels.sh
```

Prerequisites: Docker must be running and the previously built image
`msani-builder:rdkit-2025.09.1` must exist. The script does not rebuild the SDK.
Internet access is needed for pip dependencies and the clean Python images.
No GitHub Actions minutes are used.

The script builds CPython 3.10 through 3.14 wheels sequentially, repairs
them for manylinux 2.28 x86_64, and checks metadata with Twine. For each wheel
it installs into a separate `python:<version>-slim` container and runs
`python -m unittest discover -s /tmp/tests -p 'test*.py' -v`.
This includes the potentially lengthy validation-set tests.

Only tests and fixtures enter the test container, not the source package or
the SDK. Tests run on a writable temporary copy. Meeko is installed for the
PDBQT tests; Open Babel is not installed. Core native imports are checked
before Meeko is added. Runtime dependencies are resolved by pip, so their
versions can differ across Python versions; each test log records them.

Each invocation creates a new directory under `dist/linux-wheels-XXXXXXXX/`:

- `wheels/`: the repaired wheels to copy to another machine or project.
- `logs/`: separate build and unittest logs for each Python version.
- `status.txt`: pass/fail summary.
- `SHA256SUMS`: wheel checksums, relative to `wheels/`.

The script continues to the remaining Python versions after a failure and
returns a nonzero exit code if any build, metadata check, or test fails.
Failed-test wheels remain for diagnosis: do not treat their existence as a pass.
Unittest skips are reported in the logs; review them before publishing.
Nothing is uploaded to PyPI automatically.

To install on another compatible Linux x86_64 machine, select the wheel whose
`cp310`, `cp311`, `cp312`, `cp313`, or `cp314` tag matches the active Python interpreter:

```bash
python -m pip install /path/to/the/matching-wheel.whl
```

These clean-container tests do not replace testing on Ubuntu 20.04 or an
explicit matrix of supported RDKit versions.

Optional settings:

```bash
MSANI_BUILD_JOBS=4 bash packaging/linux/build_wheels.sh
```

Build a single interpreter while iterating:

```bash
MSANI_PYTHON_VERSIONS=3.12 bash packaging/linux/build_wheels.sh
```

`MSANI_BUILD_IMAGE` overrides the SDK image name. `MSANI_BUILD_RDKIT` overrides
the build-time Python RDKit version; keep it matched to the SDK in that image.
`MSANI_PYTHON_VERSIONS` selects interpreters, and `MSANI_DOCKER_PLATFORM`
selects the Docker platform (`linux/arm64` for a native aarch64 host).
Defaults are `msani-builder:rdkit-2025.09.1`, `2025.9.1`, all five versions
and `linux/amd64`.

The build and test commands themselves live in `build_in_container.sh` and
`test_in_container.sh`, and the SDK provisioning in `build_sdk.sh`, which the
Dockerfile copies into the image. GitHub Actions runs those same scripts, so a
local run exercises the CI code path.

## Build or upgrade the SDK image

Run from the repository root:

```bash
docker build --platform linux/amd64 --progress=plain \
  -t msani-builder:rdkit-2025.09.1 \
  -f packaging/linux/Dockerfile packaging/linux
bash packaging/linux/build_wheels.sh
```

Do not clear the Docker cache or use `--no-cache`. Docker rebuilds changed
layers and retains reusable layers. SDK arguments now follow the Boost layer
to avoid rebuilding Boost on future SDK upgrades. This first Dockerfile
reorganization may rebuild earlier layers once. The old image remains available.
The image tag alone does not select RDKit: the Dockerfile's `RDKIT_REF` does.
