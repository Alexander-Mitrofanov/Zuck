# my_mcaat

`my_mcaat` is an experimental, evidence-driven fork of MCAAT for improved
CRISPR array detection in raw metagenomic reads. It retains MCAAT-compatible
text output while making expensive cycle discovery replayable.

## Input modes

Exactly one input mode is required:

```bash
my_mcaat --input-files reads_R1.fastq reads_R2.fastq --output-folder result
my_mcaat --graph graph-directory --output-folder result
my_mcaat --cycles retained-cycles.txt --output-folder replay-result
```

Full runs delete the large graph and retain `cycles.txt` by default. Use
`--retain-graph` to keep the graph or `--cleanup-cycles` to remove cycle
evidence. The default maximum cycle length is 90, raised from MCAAT's 77 after
the benchmark showed a near-perfect hard miss boundary above 77 nt.

Compared with MCAAT 1.0.0, post-processing also preserves capped cycle
evidence, confines spacer collapse to bounded repeat families, prevents
single-link repeat chaining, and applies repeat-boundary correction only to
supported terminal overhangs. Numbered result files from an earlier run are
removed before new results are written, so an output directory can be reused
safely.

## Build

```bash
sudo apt-get install build-essential cmake git gettext zlib1g-dev \
  libbz2-dev libomp-dev libpthread-stubs0-dev
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

The executable is `build/my_mcaat`. The Dockerfile builds this checked-out
source rather than cloning upstream.

```bash
docker build -t my_mcaat .
docker run --rm --user "$(id -u):$(id -g)" \
  -v "$PWD/data:/data:ro" -v "$PWD/result:/out" \
  my_mcaat --input-files /data/reads_R1.fastq /data/reads_R2.fastq \
  --output-folder /out --threads 10 --ram 8G --cycle-max-length 90
```

See `DEVELOPMENT.md` for the current council decisions, benchmark gates, and
iteration history. The original MCAAT user guide is retained as `readme.md`.
