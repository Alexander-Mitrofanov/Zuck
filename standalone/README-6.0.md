# Zuck 6.0 for Linux x86-64

Zuck detects CRISPR repeats and spacers in metagenomic reads. The archive includes a portable static executable, the matching source code, a benchmark preset and a small FASTQ example.

## Run

```sh
sha256sum -c SHA256SUMS
./zuck --version
./zuck --input-files example/reads.fastq.gz --output-folder demo-result
./zuck --input-files reads_R1.fastq.gz reads_R2.fastq.gz --output-folder results
```

Supply one FASTQ for single-end reads or both mates in R1, R2 order for paired reads. Plain and gzip-compressed reads are supported. Use a fresh directory for each result.

The launcher and Docker entry point load `v6.0.conf`. This preset enables read-first discovery, adaptive graph coverage and read-supported family boundary correction, using four threads and a 4 GiB graph construction budget. Total memory use can be higher. Command-line options override the preset. For example, add `--read-first-family-recut false` to compare the spacer inventory with the longer-repeat correction switched off. Running `./bin/Zuck` directly uses native defaults.

| File in the result directory | Contents |
| --- | --- |
| `read_first_spacers.fasta` | Combined inventory of distinct spacer sequences |
| `read_first_entries.tsv` | Origin and supporting evidence for each selected spacer |
| `read_first_families.tsv` | Template and spacer-context support for each repeat family |
| `read_first_family_recut.tsv` | Evidence and decisions for candidate changes to graph boundaries |
| `supported_spacers.fasta` | Graph-supported inventory |
| `run_manifest.json` | Software version, source fingerprint, effective settings and output counts |

See [release notes](CHANGELOG.md) and the [final coverage scores](benchmark/final_spacer_performance.tsv). The [related-family controls](benchmark/recut_v6_control_scores.tsv) show where the optional boundary correction changed recovery. [Build fingerprints and validation](BUILD.txt)

## Source and Docker

`zuck-source.tar.gz` contains the release source. Use the following commands on Linux with GCC, CMake, Git, OpenMP and static zlib and bzip2 development libraries.

```sh
tar -xzf zuck-source.tar.gz
cmake -S zuck-source -B build -DCMAKE_BUILD_TYPE=Release -DSTATIC_BUILD=ON -DPORTABLE_BUILD=ON -DBUILD_TESTING=OFF
cmake --build build --target Zuck --parallel 4
```

CMake fetches pinned fmt and parallel-hashmap dependencies when they are unavailable locally. The package Dockerfile uses the prebuilt binary. Run these commands from the extracted package directory.

```sh
docker build --platform linux/amd64 -t zuck:6.0 .
docker run --rm --platform linux/amd64 --tmpfs /tmp --mount "type=bind,source=$PWD,target=/data" zuck:6.0 --input-files example/reads.fastq.gz --output-folder docker-result
```

See `LICENSE.txt`, `UPSTREAM.md` and `licenses/` for attribution and license texts.
