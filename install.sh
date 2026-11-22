#!/bin/bash

set -euo pipefail

INSTALL=false
CLEAN=false
SOURCE_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
BUILD_DIR="${SOURCE_DIR}/build"

# Parse arguments
for arg in "$@"; do
    case $arg in
        --install)
            INSTALL=true
            ;;
        -h|--help)
            echo "Usage: ./install.sh [--install] [--clean]"
            exit 0
            ;;
        --clean)
            CLEAN=true
            ;;
        *)
            echo "Unknown option: $arg"
            echo "Usage: ./install.sh [--install] [--clean]"
            exit 1
            ;;
    esac
done

if $CLEAN; then
    echo "Cleaning build directory..."
    rm -rf -- "${BUILD_DIR}"
    echo "Clean complete."
    exit 0
fi

echo "Installing system dependencies..."
sudo apt-get update
sudo apt-get install -y \
  build-essential \
  cmake \
  git \
  gettext \
  zlib1g-dev \
  libbz2-dev \
  libomp-dev \
  libpthread-stubs0-dev

echo "Creating build directory..."
mkdir -p "${BUILD_DIR}"

echo "Running CMake..."
cmake -S "${SOURCE_DIR}" -B "${BUILD_DIR}" -DCMAKE_BUILD_TYPE=Release

echo "Building..."
cmake --build "${BUILD_DIR}" -j"$(nproc)"

echo "Build complete. Binary: ${BUILD_DIR}/Zuck"

if $INSTALL; then
    echo "Installing binary to /usr/local/bin..."
    sudo cp "${BUILD_DIR}/Zuck" /usr/local/bin/
    echo "Installed. Run 'Zuck' from anywhere."
fi
