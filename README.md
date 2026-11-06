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
removed before new results are written to avoid stale files. Replacement is
not atomic, so save previous results elsewhere before reusing a directory.

## Experimental two-sided repeat boundaries

Use `--boundary-mode two-sided` to reconstruct the circular sequence before
estimating both ends of its repeat. This prevents the serialized cycle closure
from being mistaken for additional sequence, allows supported repeat variants,
and keeps each emitted repeat/spacer pair tied to an observed cycle. The
default remains `legacy`; settings files accept `boundary-mode=two-sided`.

```bash
my_mcaat --cycles retained-cycles.txt --output-folder two-sided-result \
  --boundary-mode two-sided --threads 4 --ram 2G
```

On the frozen development benchmark, strict allele-aware comparison changed
as follows. Each column compares exactly the same retained cycle input.

| Reference set | Clean arrays | Complete misses | Exact spacer matches with compatible same-row repeat |
|---|---:|---:|---:|
| CRISPRidentify Bona-fide | 15 → 22 | 20 → 7 | 655 → 1,119 |
| CRISPRCasFinder Evidence-4 | 17 → 24 | 19 → 4 | 762 → 1,568 |

The exact counts deduplicate spacers within each reference; the two callers
overlap and sometimes disagree on boundaries. All 5,384 emitted two-sided
repeat/spacer pairs reproduce an entire input circular sequence, allowing
rotation and reverse complement. Similar gains held across several retained
cycle inputs, but some small families regress and tandem-repeat backgrounds
can still produce candidates. This mode is therefore opt-in: the benchmark
does not establish biological specificity, unique loci, or spacer order.

Current graph search also fixes lost shorter paths in its reachability check,
accepts the configured maximum length, preserves 64-bit node identifiers,
and prevents repeated closure emission from using up the cycle cap. Sorted
start nodes reproduced identical cycle files in two four-thread runs. Search
can still depend on thread count and effective OpenMP team size, and the RAM
setting does not impose a hard cap on traversal memory. Fresh graph results
can change; legacy compatibility refers to replaying the same cycle input.

The current development version can additionally emit a tightly bounded,
experimental alternative-repeat view when one accepted array contains a
uniquely owned, well-supported repeat allele that its consensus cannot
represent. Enable it explicitly with `--experimental-alternative-views` (or
`experimental-alternative-views=true` in a settings file). In legacy-mode frozen-cycle
replays, primary `CRISPR_Arrays_N.txt` files match the release apart from
timestamps; supplemental views start in the next numbered file, are marked as experimental,
and contain only repeat and spacer sequences already present in their parent
array. The primary-only `getSystems()` interface excludes these views. See
`DEVELOPMENT.md` for the exact gates, benchmark result, and current
held-out-validation caveat. The view repeats existing evidence; a later audit
found that the original benchmark had missed the parent allele during
comparison, so its score change does not establish discovery of a new spacer.

Cycle evidence is stored in a lossless compact overlap representation. On the
frozen benchmark this reduced `cycles.txt` from about 422 MB to 23 MB while
replaying to the same primary arrays. Explicit RAM values are interpreted as
GiB and converted to bytes for graph construction; automatic memory selection
uses the process' available host/cgroup budget.

## Build

```bash
sudo apt-get install build-essential cmake git gettext zlib1g-dev \
  libbz2-dev libomp-dev libpthread-stubs0-dev
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

The executable is `build/my_mcaat`. The Dockerfile builds this checked-out
source rather than cloning upstream.

To build and run the tests, configure with `-DBUILD_TESTING=ON`, build the
project, then run `ctest --test-dir build --output-on-failure`. The CLI I/O
and boundary-mode tests, plus the production-enumerator oracle test, are
registered when Python 3 is available.

```bash
docker build -t my_mcaat .
docker run --rm --user "$(id -u):$(id -g)" \
  -v "$PWD/data:/data:ro" -v "$PWD/result:/out" \
  my_mcaat --input-files /data/reads_R1.fastq /data/reads_R2.fastq \
  --output-folder /out --threads 10 --ram 8G --cycle-max-length 90
```

See `DEVELOPMENT.md` for the current council decisions, benchmark gates, and
iteration history. The original MCAAT user guide is retained as `readme.md`.
