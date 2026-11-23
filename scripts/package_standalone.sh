#!/usr/bin/env bash
# Assemble the two built executables and their runtime files in a new directory.
set -euo pipefail
src=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
build=${1:?supply the CMake build directory}
dest=${2:?supply a new bundle directory}
mkdir -- "$dest"
install -m 0755 "$build/Zuck" "$build/zuck-kmer-spectrum" "$dest/"
install -m 0644 "$src"/{LICENSE.txt,UPSTREAM.md} "$dest/"
install -m 0644 "$src"/standalone/{README.md,ADAPTIVE.md,COVERAGE.md,v7.0.conf,v6.1.conf,zuck_depth.py,zuck_depth_model.py,zuck_sparse_assisted.py} "$dest/"
(cd -- "$dest" && sha256sum -- * > SHA256SUMS)
