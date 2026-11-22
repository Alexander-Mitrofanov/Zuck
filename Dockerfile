# Stage 1: build a portable, dependency-free executable and its tests.
FROM debian:bookworm AS builder

ENV DEBIAN_FRONTEND=noninteractive
ARG BUILD_JOBS=4

RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    git \
    gettext \
    zlib1g-dev \
    libbz2-dev \
    libomp-dev \
    libpthread-stubs0-dev \
    && rm -rf /var/lib/apt/lists/*

# Python enables the CLI, result-bundle, boundary, and standalone CTest cases.
RUN apt-get update && apt-get install -y python3 time && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /app

# Build the checked-out Zuck source, including its vendored dependencies. The
# static and portable options make this exact binary suitable for both the
# scratch container and direct execution on another Linux distribution.
COPY . .

RUN cmake -S . -B build \
      -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_EXPORT_COMPILE_COMMANDS=ON \
      -DBUILD_TESTING=ON \
      -DSTATIC_BUILD=ON \
      -DPORTABLE_BUILD=ON && \
    cmake --build build --parallel "${BUILD_JOBS}"

# Stage 2: run the complete native/CLI suite against the artifact that ships.
FROM builder AS tested

ARG TEST_JOBS=4

RUN ctest --test-dir build --output-on-failure --parallel "${TEST_JOBS}" && \
    ! grep -Fq -- '-march=native' build/compile_commands.json && \
    ! readelf -l build/Zuck | grep -q 'Requesting program interpreter' && \
    ! readelf -d build/Zuck | grep -q '(NEEDED)' && \
    mkdir -p /standalone && \
    install -m 0755 build/Zuck /standalone/Zuck && \
    install -m 0644 LICENSE.txt /standalone/LICENSE.txt && \
    install -m 0644 UPSTREAM.md /standalone/UPSTREAM.md && \
    install -m 0644 standalone/README.md /standalone/README.md && \
    cd /standalone && \
    sha256sum Zuck LICENSE.txt UPSTREAM.md README.md > SHA256SUMS

# Export with:
# docker build --target standalone --output type=local,dest=dist/zuck-standalone .
FROM scratch AS standalone

COPY --from=tested /standalone/ /

ENTRYPOINT ["/Zuck"]

# The normal image contains the same tested standalone executable. A scratch
# runtime prevents accidental dependence on files that are absent from the
# exported bundle.
FROM scratch AS runtime

COPY --from=tested /standalone/Zuck /usr/local/bin/Zuck

ENTRYPOINT ["/usr/local/bin/Zuck"]
