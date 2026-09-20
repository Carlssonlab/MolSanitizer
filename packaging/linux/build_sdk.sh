#!/usr/bin/env bash
# Provision the pinned RDKit C++ SDK inside a manylinux_2_28 environment.
#
# Single source of truth for the SDK: the Dockerfile invokes the stages
# separately so that a Boost layer survives an RDKit bump, while CI calls
# "all" directly inside the upstream manylinux container. Works on x86_64
# and aarch64; nothing here is architecture specific.
set -euo pipefail

RDKIT_REF=${RDKIT_REF:-Release_2025_09_1}
RDKIT_PREFIX=${RDKIT_PREFIX:-/opt/msani-rdkit}
BOOST_PREFIX=${BOOST_PREFIX:-/opt/boost}
BOOST_VERSION=${BOOST_VERSION:-1.85.0}
BOOST_SHA256=${BOOST_SHA256:-be0d91732d5b0cc6fbb275c7939974457e79b54d6f07ce2e3dfdd68bef883b0b}
JOBS=${MSANI_BUILD_JOBS:-$(nproc)}

stage_deps() {
    dnf -y install git cmake ninja-build eigen3-devel \
        curl tar gzip zlib-devel bzip2-devel xz-devel libzstd-devel
    dnf clean all
}

# Pinned, position-independent static Boost. RDKit and MolSanitizer must agree
# on one Boost; the manylinux image does not ship a suitable development copy.
stage_boost() {
    local archive=/tmp/boost.tar.gz
    local source="/tmp/boost_${BOOST_VERSION//./_}"
    curl -fL --retry 3 \
        "https://archives.boost.io/release/${BOOST_VERSION}/source/boost_${BOOST_VERSION//./_}.tar.gz" \
        -o "$archive"
    echo "${BOOST_SHA256}  ${archive}" | sha256sum -c -
    tar -xzf "$archive" -C /tmp
    cd "$source"
    ./bootstrap.sh --prefix="$BOOST_PREFIX" \
        --with-libraries=serialization,iostreams,filesystem,system,thread,regex
    ./b2 -j"$JOBS" install \
        variant=release link=static threading=multi \
        cxxflags="-fPIC" \
        --layout=system
    cd /
    rm -rf "$archive" "$source"
}

# Static RDKit avoids a runtime dependency on an RDKit C++ shared-library ABI.
# Python wrappers are not needed: pip installs the separately pinned Python API.
stage_rdkit() {
    rm -rf /tmp/rdkit /tmp/rdkit-build
    git clone --depth 1 --branch "$RDKIT_REF" \
        https://github.com/rdkit/rdkit.git /tmp/rdkit
    cmake -S /tmp/rdkit -B /tmp/rdkit-build -G Ninja \
        -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_INSTALL_PREFIX="$RDKIT_PREFIX" \
        -DCMAKE_PREFIX_PATH="$BOOST_PREFIX" \
        -DRDK_INSTALL_INTREE=OFF \
        -DRDK_BUILD_STATIC_LIBS_ONLY=ON \
        -DRDK_BUILD_PYTHON_WRAPPERS=OFF \
        -DRDK_BUILD_CPP_TESTS=OFF \
        -DRDK_BUILD_INCHI_SUPPORT=OFF \
        -DCMAKE_POSITION_INDEPENDENT_CODE=ON \
        -DRDK_BUILD_FREETYPE_SUPPORT=OFF \
        -DRDK_INSTALL_COMIC_FONTS=OFF \
        -DRDK_BUILD_MAEPARSER_SUPPORT=OFF \
        -DRDK_BUILD_PUBCHEMSHAPE_SUPPORT=OFF
    cmake --build /tmp/rdkit-build --parallel "$JOBS"
    cmake --install /tmp/rdkit-build
    # Record provenance so a restored SDK tarball can be matched to its pin.
    install -d "$RDKIT_PREFIX/share/msani-sdk"
    git -C /tmp/rdkit rev-parse HEAD > "$RDKIT_PREFIX/share/msani-sdk/rdkit-commit.txt"
    install -m 0644 /tmp/rdkit/license.txt "$RDKIT_PREFIX/share/msani-sdk/RDKit-license.txt"
    printf '{"rdkit_ref":"%s","boost_version":"%s","arch":"%s","linkage":"static"}\n' \
        "$RDKIT_REF" "$BOOST_VERSION" "$(uname -m)" \
        > "$RDKIT_PREFIX/share/msani-sdk/msani-sdk.json"
    rm -rf /tmp/rdkit /tmp/rdkit-build
}

case "${1:-all}" in
    deps) stage_deps ;;
    boost) stage_boost ;;
    rdkit) stage_rdkit ;;
    all) stage_deps; stage_boost; stage_rdkit ;;
    *) printf 'Usage: %s [deps|boost|rdkit|all]\n' "$0" >&2; exit 2 ;;
esac
