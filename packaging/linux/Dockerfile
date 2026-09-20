# Build the private RDKit C++ SDK before building the MolSanitizer wheel.
# The result is an intermediate build artifact, never an uploaded package.
#
# The provisioning logic lives in build_sdk.sh so that this image and the
# GitHub Actions SDK job build byte-identical prefixes. Stages are invoked
# separately here to keep the expensive Boost layer cached across RDKit bumps.
ARG MANYLINUX_IMAGE=quay.io/pypa/manylinux_2_28_x86_64
FROM ${MANYLINUX_IMAGE}

COPY build_sdk.sh /opt/msani/build_sdk.sh

ARG MSANI_BUILD_JOBS=2
ENV MSANI_BUILD_JOBS=${MSANI_BUILD_JOBS}

RUN bash /opt/msani/build_sdk.sh deps

ARG BOOST_VERSION=1.85.0
ARG BOOST_PREFIX=/opt/boost
ENV BOOST_VERSION=${BOOST_VERSION} BOOST_PREFIX=${BOOST_PREFIX}
RUN bash /opt/msani/build_sdk.sh boost

# Keep SDK arguments below Boost so future SDK changes preserve its cache.
ARG RDKIT_REF=Release_2025_09_1
ARG RDKIT_PREFIX=/opt/msani-rdkit
ENV RDKIT_REF=${RDKIT_REF} RDKIT_PREFIX=${RDKIT_PREFIX}
RUN bash /opt/msani/build_sdk.sh rdkit

# Environment search paths use ':' on Linux (CMake -D lists use ';').
ENV CMAKE_PREFIX_PATH="${RDKIT_PREFIX}:${BOOST_PREFIX}"
ENV MSANI_RDKIT_LINKAGE=static
