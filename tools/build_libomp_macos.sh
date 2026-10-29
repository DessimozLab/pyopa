#!/usr/bin/env bash
# Builds LLVM's OpenMP runtime (libomp) from source and installs it into the
# given prefix. Used for the macOS wheels: unlike Homebrew's libomp, which
# targets the macOS version of the build machine, this honours
# MACOSX_DEPLOYMENT_TARGET (set in pyproject.toml), so the bundled library
# does not raise the minimum macOS version of the wheel.
set -euo pipefail

PREFIX=${1:?usage: $0 <install-prefix>}
LLVM_VERSION=23.1.4
LLVM_SHA256=819ecf7351f1a66f76b2d701d8942faaf5faeb324cb90c9d5b090de142434bb2
TARGET=${MACOSX_DEPLOYMENT_TARGET:-11.0}
echo "Building libomp ${LLVM_VERSION} for macOS ${TARGET} ($(uname -m))"

BUILD_DIR=$(mktemp -d)
trap 'rm -rf "$BUILD_DIR"' EXIT
cd "$BUILD_DIR"

TARBALL=llvm-project-${LLVM_VERSION}.src.tar.xz
curl -fsSL -o "$TARBALL" \
    "https://github.com/llvm/llvm-project/releases/download/llvmorg-${LLVM_VERSION}/${TARBALL}"
echo "${LLVM_SHA256}  ${TARBALL}" | shasum -a 256 -c -
# only the parts needed for a runtimes build of openmp
SRC=llvm-project-${LLVM_VERSION}.src
tar xf "$TARBALL" "$SRC/openmp" "$SRC/runtimes" "$SRC/cmake" "$SRC/llvm/cmake" \
    "$SRC/llvm/utils/lit" "$SRC/llvm/utils/llvm-lit" "$SRC/third-party/unittest"

cmake -S "$SRC/runtimes" -B build \
    -DLLVM_ENABLE_RUNTIMES=openmp \
    -DLLVM_INCLUDE_TESTS=OFF \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX="$PREFIX" \
    -DCMAKE_OSX_ARCHITECTURES="$(uname -m)" \
    -DCMAKE_OSX_DEPLOYMENT_TARGET="$TARGET" \
    -DCMAKE_INSTALL_NAME_DIR="$PREFIX/lib" \
    -DLIBOMP_INSTALL_ALIASES=OFF \
    -DOPENMP_ENABLE_LIBOMPTARGET=OFF \
    -DOPENMP_ENABLE_OMPT_TOOLS=OFF \
    -DLIBOMP_OMPD_SUPPORT=OFF
cmake --build build --parallel "$(sysctl -n hw.ncpu)"
cmake --install build
