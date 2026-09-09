# Build the private RDKit C++ SDK before building the MolSanitizer wheel.
# The result is an intermediate build artifact, never an uploaded package.
FROM quay.io/pypa/manylinux_2_28_x86_64

RUN dnf -y install git cmake ninja-build eigen3-devel \
        curl tar gzip zlib-devel bzip2-devel xz-devel libzstd-devel \
    && dnf clean all

# Build a pinned, modern Boost with position-independent static libraries.
RUN curl -fL --retry 3 \
        https://archives.boost.io/release/1.85.0/source/boost_1_85_0.tar.gz \
        -o /tmp/boost.tar.gz \
    && echo "be0d91732d5b0cc6fbb275c7939974457e79b54d6f07ce2e3dfdd68bef883b0b  /tmp/boost.tar.gz" \
        | sha256sum -c - \
    && tar -xzf /tmp/boost.tar.gz -C /tmp \
    && cd /tmp/boost_1_85_0 \
    && ./bootstrap.sh --prefix=/opt/boost \
        --with-libraries=serialization,iostreams,filesystem,system,thread,regex \
    && ./b2 -j2 install \
        variant=release link=static threading=multi \
        cxxflags="-fPIC" \
        --layout=system

# Keep SDK arguments below Boost so future SDK changes preserve its cache.
ARG RDKIT_REF=Release_2025_09_1
ARG RDKIT_PREFIX=/opt/msani-rdkit

# Environment search paths use ':' on Linux (CMake -D lists use ';').
ENV CMAKE_PREFIX_PATH="${RDKIT_PREFIX}:/opt/boost"

RUN git clone --depth 1 --branch "${RDKIT_REF}" \
        https://github.com/rdkit/rdkit.git /tmp/rdkit

# Static RDKit avoids a runtime dependency on an RDKit C++ shared-library ABI.
# Python wrappers are not needed: pip installs the separately pinned Python API.
RUN cmake -S /tmp/rdkit -B /tmp/rdkit-build -G Ninja \
        -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_INSTALL_PREFIX="${RDKIT_PREFIX}" \
        -DCMAKE_PREFIX_PATH=/opt/boost \
        -DRDK_INSTALL_INTREE=OFF \
        -DRDK_BUILD_STATIC_LIBS_ONLY=ON \
        -DRDK_BUILD_PYTHON_WRAPPERS=OFF \
        -DRDK_BUILD_CPP_TESTS=OFF \
        -DRDK_BUILD_INCHI_SUPPORT=OFF \
        -DCMAKE_POSITION_INDEPENDENT_CODE=ON \
        -DRDK_BUILD_FREETYPE_SUPPORT=OFF \
        -DRDK_INSTALL_COMIC_FONTS=OFF \
        -DRDK_BUILD_MAEPARSER_SUPPORT=OFF \
        -DRDK_BUILD_PUBCHEMSHAPE_SUPPORT=OFF \
    && cmake --build /tmp/rdkit-build --parallel 2 \
    && cmake --install /tmp/rdkit-build

ENV MSANI_RDKIT_LINKAGE=static
