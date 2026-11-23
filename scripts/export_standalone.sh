#!/usr/bin/env bash
set -euo pipefail

source_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
output=${1:-"${source_dir}/dist/zuck-standalone"}
build_jobs=${BUILD_JOBS:-4}
test_jobs=${TEST_JOBS:-4}
image=zuck-standalone-export:local
container=

fail() {
    echo "Error: $*" >&2
    exit 1
}

cleanup() {
    if [[ -n "${container}" ]]; then
        docker rm -f "${container}" >/dev/null 2>&1 || true
    fi
}
trap cleanup EXIT
trap 'exit 130' HUP INT TERM

[[ "${build_jobs}" =~ ^[1-9][0-9]*$ ]] || fail "BUILD_JOBS must be a positive integer"
[[ "${test_jobs}" =~ ^[1-9][0-9]*$ ]] || fail "TEST_JOBS must be a positive integer"
[[ ! -e "${output}" && ! -L "${output}" ]] || fail "output already exists: ${output}"
command -v docker >/dev/null || fail "docker is required"

docker build \
    --target standalone \
    --build-arg "BUILD_JOBS=${build_jobs}" \
    --build-arg "TEST_JOBS=${test_jobs}" \
    -t "${image}" \
    "${source_dir}"

container=$(docker create "${image}")
mkdir -p -- "${output}"
docker cp "${container}:/SHA256SUMS" "${output}/SHA256SUMS"
while read -r digest artifact; do
    docker cp "${container}:/${artifact}" "${output}/${artifact}"
done < "${output}/SHA256SUMS"

(
    cd -- "${output}"
    sha256sum -c SHA256SUMS
    ./Zuck --help >/dev/null
)

echo "Standalone bundle: ${output}"
