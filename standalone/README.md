# Zuck standalone for Linux

`Zuck` is the exact executable tested and installed in the accompanying Docker
image. It is statically linked, so it does not require Docker, a C++ runtime,
OpenMP, zlib, or other userspace libraries on the destination machine.

The executable is portable across Linux distributions on the CPU architecture
for which the image was built. It deliberately avoids `-march=native`; the
kernel and CPU must still match the build architecture (for example,
`linux/amd64`).

Verify and run it from the extracted directory:

```bash
sha256sum -c SHA256SUMS
./Zuck --help
./Zuck --input-files reads_R1.fastq reads_R2.fastq \
  --output-folder result --threads 4 --ram 4G
```

To build a tested bundle from the source tree with Docker BuildKit:

```bash
scripts/export_standalone.sh
```

The Docker build is capped at four compilation and test workers by default.
Override `BUILD_JOBS` or `TEST_JOBS` only when the build environment permits it.

The optional coverage preflight estimates sequencing depth from FASTQ and can select a Zuck 6.1 calling mode. [Build and run the preflight](COVERAGE.md)
