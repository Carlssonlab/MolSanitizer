# Binary-wheel packaging

`linux/Dockerfile` builds the pinned RDKit 2025.09.5 C++ SDK with static
archives and installs it at `/opt/msani-rdkit`. Build MolSanitizer inside that
image with `-DMSANI_RDKIT_LINKAGE=static`; do not use a developer Conda prefix
for an uploadable wheel.

Before publishing, run `auditwheel show` and `auditwheel repair` on the raw
wheel, then install the repaired wheel into a clean `python:3.12-slim` image
with only pip-installed runtime dependencies. The test must import all three
native modules and generate a conformer. It must not mount or reference a
Conda installation.

This is the first Linux x86_64 / CPython 3.12 packaging target. Each Python
and platform wheel needs its own build and clean-install test.

For the current local prototype, `/home/phonglam/msani-rdkit-sdk` is a valid
RDKit 2025.09.5 SDK and the build selects RDKit's `*_static` imported targets.
The resulting WSL-built wheel is limited to `manylinux_2_39_x86_64`; use the
manylinux container for the wider `manylinux_2_28_x86_64` production target.

## Local platform wheel matrices

- Linux x86_64: `linux/build_wheels.sh` (Docker).
- Windows x64: [Windows instructions](windows/README.md).
- macOS arm64 and x86_64: [macOS instructions](macos/README.md), or run
  `bash packaging/macos/build_wheels.sh` on a native Mac of that architecture.

## Shared build steps

The SDK provisioning and the per-interpreter build and test commands live in
scripts that both local runs and GitHub Actions invoke, so CI wheels and local
wheels come from identical commands:

| Script | Used by |
| --- | --- |
| `linux/build_sdk.sh` | `linux/Dockerfile` (per stage) and the SDK workflow |
| `linux/build_in_container.sh` | `linux/build_wheels.sh` and the wheel workflow |
| `linux/test_in_container.sh` | the same two, in a clean `python:X.Y-slim` image |
| `macos/build_sdk.sh` | `macos/build_wheels.sh` and the SDK workflow |
| `fetch_sdk.sh`, `upload_sdk.sh` | the workflows, to move SDK release assets |

## Continuous integration

See [RELEASE.md](RELEASE.md) for the two workflows and the manual upload step.
Wheels are never published automatically.
