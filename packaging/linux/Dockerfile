# Build the private RDKit C++ SDK before building the MolSanitizer wheel.
# The result is an intermediate build artifact, never an uploaded package.
FROM quay.io/pypa/manylinux_2_28_x86_64

ARG RDKIT_REF=Release_2025_09_5
ARG RDKIT_PREFIX=/opt/msani-rdkit

RUN dnf -y install git cmake ninja-build boost-devel eigen3-devel && dnf clean all

RUN git clone --depth 1 --branch "${RDKIT_REF}" \
        https://github.com/rdkit/rdkit.git /tmp/rdkit

# Static RDKit avoids a runtime dependency on an RDKit C++ shared-library ABI.
# Python wrappers are not needed: pip installs the separately pinned Python API.
RUN cmake -S /tmp/rdkit -B /tmp/rdkit-build -G Ninja \
        -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_INSTALL_PREFIX="${RDKIT_PREFIX}" \
        -DBUILD_SHARED_LIBS=OFF \
        -DRDK_BUILD_PYTHON_WRAPPERS=OFF \
        -DRDK_BUILD_CPP_TESTS=OFF \
        -DRDK_BUILD_INCHI_SUPPORT=OFF \
    && cmake --build /tmp/rdkit-build --parallel \
    && cmake --install /tmp/rdkit-build

ENV CMAKE_PREFIX_PATH=${RDKIT_PREFIX}
ENV MSANI_RDKIT_LINKAGE=static
