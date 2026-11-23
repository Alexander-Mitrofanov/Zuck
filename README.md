# Zuck 7.0

Zuck detects CRISPR repeat families and spacer sequences in metagenomic reads. Version 7.0 estimates sequencing depth directly from FASTQ. When reads support a uniform estimate below 5×, it searches for complete spacers observed only once whose boundaries are supported by repeats on both sides. At intermediate, high, mixed or uncertain coverage, it uses the full 6.1 calling logic.

Download the Linux x86-64 bundle, extract it and run the following commands from its directory. The adaptive launcher requires Python 3. [Zuck 7.0 release](https://github.com/Alexander-Mitrofanov/Zuck/releases/tag/v7.0)

```sh
sha256sum -c SHA256SUMS
./Zuck --version
python3 zuck_depth.py adaptive --input-files R1.fastq.gz R2.fastq.gz --output-folder result
```

Read `result/read_first_spacers.fasta` for the native spacer inventory and `result/coverage_routing.json` for the depth estimate and calling mode. Native Zuck runs independently of CRASS. The optional assisted mode uses an external CRASS installation and writes its own inventory. [Calling modes and outputs](standalone/ADAPTIVE.md) · [Exact and tolerant spacer F1](docs/benchmark-7.0.md)

## Build and test

```bash
sudo apt-get install build-essential cmake git gettext zlib1g-dev \
  libbz2-dev libomp-dev libpthread-stubs0-dev
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON
cmake --build build -j
ctest --test-dir build --output-on-failure
```

The native caller is `build/Zuck` and the coverage scanner is `build/zuck-kmer-spectrum`. CMake fetches the pinned parallel-hashmap and
fmt dependencies when they are not already available. Python 3.11 or newer enables the
CLI, result-bundle, boundary, and cycle-enumerator tests. The Dockerfile builds
the checked-out source and its vendored dependencies.

### Standalone Linux bundle

The Docker build compiles the portable static caller and scanner. The default scratch image runs the native caller directly. To export the adaptive bundle for use with Python on the host, run the following commands.

```bash
scripts/export_standalone.sh
(cd dist/zuck-standalone && sha256sum -c SHA256SUMS)
dist/zuck-standalone/Zuck --help
```

The export script works with both the current BuildKit builder and Docker's
legacy builder. With BuildKit, the equivalent direct export is
`docker build --target standalone --output
type=local,dest=dist/zuck-standalone .`.

The bundle includes both executables, the Python launcher and model, presets, checksums, license and upstream attribution. Both executables are statically linked and built without host-specific `-march=native` instructions.
Its CPU architecture still matches the Docker build platform; use an explicit
`--platform` when producing architecture-specific releases.

## Run

Exactly one input mode is required:

```bash
build/Zuck --input-files reads_R1.fastq reads_R2.fastq \
  --output-folder result --threads 8 --ram 4G
build/Zuck --graph graph-directory --output-folder graph-result
build/Zuck --cycles result/cycles.txt --output-folder replay-result \
  --grouping-reads reads_R1.fastq reads_R2.fastq --threads 8 --ram 4G
```

Full runs retain `cycles.txt` and remove the graph by default. Use
`--retain-graph` to keep the graph or `--cleanup-cycles` to remove retained
cycles. Supplied graphs are preserved. Cleanup removes only recognized generated
files from graph folders and leaves unrelated files intact. The maximum cycle length defaults to 90 bases. Explicit RAM values
use binary units; `4G` means 4 GiB. The RAM setting budgets graph construction
but does not impose a hard cap on traversal memory. Automatic memory
selection uses the available host/cgroup budget.

Settings files accept `key=value` entries corresponding to CLI options, such
as `boundary-mode=two-sided`, `threads=8`, and `ram=4G`. CLI values override
file settings. Numeric and boolean values are validated in both interfaces;
malformed values fail before output is changed. Quote individual filenames
containing spaces in `input-files` or `grouping-reads`, for example
`input-files="reads/read one.fastq" "reads/read two.fastq"`. Scalar path settings
such as `output-folder` take the rest of the value literally, without quotes. Run `build/Zuck --help` for the complete interface.

### Read-first spacer discovery

Use `--read-first true` with single-end or paired FASTQ reads and an enabled
call profile. Zuck discovers candidate repeats directly in individual reads,
pools evidence across spacers that share a repeat, and checks both
repeat–spacer junctions. Once other spacers establish the family, a single
well-bounded observation can support a complete spacer.

Use `--adaptive-coverage true` to let graph-edge multiplicities guide cycle
search. When repeated edges have low support, Zuck adjusts the seed cutoff
and retains singleton edges. Explicit settings for either parameter take
precedence. The `coverage_support` section of `run_manifest.json` records the
observed multiplicities and the decision.

The experimental `--read-first-family-recut true` option requires read-first
discovery. When two distinct repeat families share a core, a longer repeat
observed in reads can correct a shorter graph boundary. For an ambiguous
placement, the direct check requires at least two distinct fragments carrying
the exact repeat–spacer–repeat sequence at Q20 or higher. Other fragments
must independently establish the competing repeat families.

The final coverage benchmark used the following native command with paired reads.

```bash
build/Zuck --input-files reads_R1.fastq reads_R2.fastq --output-folder result \
  --threads 4 --ram 4G --call-profile standard --read-recruitment true \
  --local-boundary-graph false --motif-cliff-boundary false \
  --singleton-precision-prune true --circular-unit-dedup false \
  --directed-primary-groups true --cross-variant-scaffolds true \
  --adaptive-coverage true --read-first true --read-first-family-recut true
```

| File in the output folder | What it contains |
| --- | --- |
| `read_first_spacers.fasta` | Combined graph and read-first inventory of distinct, strand-canonical spacers. This is the FASTA used for the benchmark. |
| `read_first_entries.tsv` | The selected spacer, its repeat family, and the graph or read evidence that supplied it. |
| `read_first_families.tsv` | Counts of supporting templates and distinct spacer contexts for each repeat. A repeat family can occur at several loci. |
| `read_first_family_recut.tsv` | Each boundary decision and its competing repeat when family recutting is enabled. `full_pairs` counts unambiguous full-unit templates in the extractor. `raw_q20_full_templates` counts distinct fragments carrying the exact full unit in the direct check, including reads compatible with a shorter boundary. |

#### Coverage benchmark

The frozen 6.0 release candidate and CRASS 1.0.1 were compared on 30
paired-end InSilicoSeq HiSeq libraries from three held-out five-genome
communities. At each coverage, F1 was pooled across both draws of all three
communities against the union of CRISPRidentify Bona-fide and CRISPRCasFinder
Evidence 4 spacers. Tolerant matching required at least 90% identity and
90% reciprocal span. Mixed libraries assigned different genomes depths
ranging from 1× to 20×.

| Coverage | Exact F1, Zuck | Exact F1, CRASS | Tolerant F1, Zuck | Tolerant F1, CRASS |
| --- | ---: | ---: | ---: | ---: |
| 1× | 39.82% | 36.52% | 41.97% | 44.38% |
| 3× | 69.71% | 64.51% | 71.51% | 74.41% |
| 8× | 80.80% | 67.83% | 83.04% | 84.63% |
| 20× | 79.78% | 70.44% | 82.17% | 84.78% |
| Mixed | 67.79% | 55.92% | 70.20% | 71.46% |

Zuck led at all five depths on exact F1 and CRASS led at all five on tolerant
F1. At 20×, Zuck recovered 1,341 tolerant matches from 1,730 predictions.
CRASS recovered 1,259 from 1,436. Turning the family boundary correction
off gave identical inventories in these 30 libraries. On a separate targeted
80× pool of three loci with related repeats, enabling it increased exact
recovery from 46 to 199 of 223 annotated spacers. CRASS recovered 217.

In the preceding cohort, median template support among detected, uniquely
assignable repeat types rose from 6 at 1× to 31 at 20×. Repeat copy number,
spacer diversity, read geometry and sampling limits prevent direct conversion
of these counts to absolute genome coverage.

### Experimental read recruitment

Add `--read-recruitment true` to a full read run or a cycle replay with
`--grouping-reads` to recruit complete spacers around existing graph repeat
seeds. The pass supports adjacent repeat mutations and clipped repeat
junctions while preserving exact spacer alleles and original row IDs. It is
disabled by default and requires two-sided or branch-aware boundaries.

The additional evidence channel has explicit provenance, a physical-triplet
screen for tandem-repeat explanations, and a conservative quality/depth screen
for one-substitution allele ambiguity. Ambiguous observed
alleles remain inclusive. Automatic recruitment promotion requires an original
graph row and read-supported diverse spacer context in every owning family.
Newly observed pairs remain unresolved, including genuine novel alleles.
This pass does not recover seedless families or resolve an already ambiguous
repeat boundary. See [result interpretation](docs/results.md#graph-seeded-read-recruitment)
for bounds, selection rules, evidence tables and statistical limitations.

### Search reproducibility and optional length limits

Cycle traversal uses eight logical search partitions by default. Their sorted,
contiguous start-node ranges and visited histories are independent of the
requested or effective OpenMP worker count. `--threads` controls execution;
`--traversal-partitions` controls the search histories and can change the
discovered evidence. Eight partitions preserve the historical eight-worker
development configuration. Enumeration uses at most one worker per logical
partition (eight with the default); graph building and start-node scanning
can still use the full requested team. To reproduce an earlier run using another worker
count, set the partition count to that earlier effective team size. This
does not make the bounded, visited-history-dependent search exhaustive.

`--low-abundance true` now retains graph edges with multiplicity one. The
default is `false`, preserving the previous effective singleton removal:
earlier versions removed these edges even when the flag was true. Tip pruning,
start-node multiplicity thresholds and other search limits still apply.
An old settings file explicitly containing `low-abundance=true` therefore
changes behavior after this fix; use `false` to preserve its old effect.

`--spacer-min-length`, `--spacer-max-length`, `--repeat-min-length`, and
`--repeat-max-length` are optional, inclusive **output constraints**, each
defaulting to zero (no additional bound). The same settings-file keys now
take effect; their previous nominal defaults were unused. Active constraints
filter final de novo and database rows and group repeat consensuses before
assigning row IDs and producing linkage/confidence sidecars. A filtered group
must retain at least two distinct exact spacers. These options do not expand
discovery limits or rerun boundary inference; for example, lowering a spacer
minimum cannot recover candidates excluded by the internal 20-base floor.

## Boundary modes

| Mode | Behavior |
|---|---|
| `legacy` | Historical boundary inference, available for reproducibility. |
| `two-sided` | Reconstructs circular sequence, infers both repeat ends, and preserves observed repeat/spacer pairing. |
| `branch-aware` | Adds experimental separation of divergent repeat families sharing an anchor. |

The default is `two-sided`, selected for the strongest measured recovery and
concordance in the completed 12 November development benchmark. Explicit
`--boundary-mode legacy` remains available for historical replay. Branch-aware
mode improves several synthetic mixtures but does not improve high-confidence
reference recovery in that benchmark.

`--boundary-mode branch-aware --repeat-mixture-refinement true` additionally
enables an experimental model for nested repeat-family mixtures. It permits
five-base signatures on one repeat end, up to two partition levels, and at
least three edits between inferred branch repeats. Each branch still needs
three distinct circular cores; at most four branches must cover at least 80%
of cores, and a replacement must extend the repeat by at least five bases.
The default is **false**, and all historical modes remain unchanged when off.
Settings files accept `repeat-mixture-refinement`; the run manifest and text
output identify the experimental model.

This rule was developed after inspecting the 15 November benchmark and is not
independent validation. It repairs tested repeat mixtures, but **can assign
conserved spacer termini to the repeat**: two adversarial controls containing
40 true spacers each fell from 40 exact calls to zero, while retaining 40
read-supported shifted calls. Exact raw-read support cannot determine an
otherwise ambiguous repeat/spacer boundary. Do not enable this option as a
general accuracy upgrade; retain baseline output when evaluating it.

Modern boundary modes retain small, equally supported terminal spacer variants
by default (`--terminal-variant-rescue false` disables this). A candidate must
have the same traversal support as the selected representative, differ only
by one terminal base, and belong to a component of at most eight observed
alleles with representative support of at least three. Only observed
repeat/spacer pairs are added; all previously accepted rows remain present.
Traversal counts are not independent read support.

Read-supported grouping is enabled by default (`--read-linked-groups false`
disables it). It adds supplementary content subsets using exact repeat and
spacer-flank junctions in the input reads. Every accepted edge requires two
distinct full sequence units and two junction-start offsets; identical read
pairs cannot supply both units. These counts do not establish independent
molecules. Across each exact repeat/flank span, all bases must have Phred
quality at least 10 and the sum of their error probabilities must be at most
one. Only complete unbranched components with at least three spacer
types are emitted, and every parent group and oriented row is retained.
FASTQ qualities determine support; FASTA records cannot provide this support.
Cycles and graph replays need `--grouping-reads` to reproduce these views.
Legacy mode does not apply terminal rescue or read-supported grouping.

The selected default recovers 1,121/1,517 exact CI spacers and 1,570/1,703 CCF
spacers, two more per panel than the prior default. Read-linked views improve
strict clean-array recovery from 25 to 36 of 58 CI arrays and from 27 to 36 of
59 CCF arrays. All 5,421 original oriented rows remain, with nine added rows
in the same 610 primary groups. Every previously source-supported spacer and
literal paired sequence is retained. See `DEVELOPMENT.md` for full-run
measurements, rejected alternatives, and limits.

The optional `--short-spacer-filter true` setting rejects groups explained by a single
observed longest spacer plus shorter prefix/suffix variants. Every shorter
variant must be below 23 bases and differ by at most two terminal bases from
that anchor. These fixed thresholds address the 23-mer collapse behavior;
`spacer-min-length` does not change them. The filter applies to two-sided and
branch-aware modes; legacy mode does not apply it. The default is **false**;
settings files can enable it with `short-spacer-filter=true`. In the
development ablation it removed eight distinct spacers, all exactly supported
by the source genomes. Its F1 gain against incomplete reference panels does
not justify making that loss the default. Removing such groups does not
establish that they were false biological calls.

The optional `--experimental-alternative-views` flag emits bounded alternate
repeat views of already accepted evidence. It does not discover new spacers.
These views are labelled separately and are excluded from the primary-only
`getSystems()` interface.

## Candidate evidence before spacer collapse

`--candidate-evidence true` writes `candidate_alleles.tsv` for every exact
oriented repeat/spacer pair extracted by two-sided or branch-aware boundary
inference, before repeat-family and spacer collapse. This optional diagnostic
adds a read scan and does not change accepted calls or their confidence.
Replay requires `--grouping-reads` for read evidence; without reads the
corresponding scores are `NA`.

The table records traversal observations, distinct canonical circular cores,
quality-qualified repeat/spacer/context read scores, bounded context sequence
units and start offsets, saturation and observability. It also marks whether
the exact pair was emitted in the primary output and whether it passes the
explicit pair-length constraints. That length flag does not test group-level
consensus or minimum-content requirements. Candidate IDs follow sorted exact
pairs and are stable for a fixed extracted candidate set.

Candidates are not additional CRISPR calls. They can include incorrect cuts,
sequencing errors and redundant representations; even read-supported context
does not establish a biological repeat/spacer boundary. Traversal multiplicity
is not read depth. The table excludes database-only candidates and candidates
already lost during graph traversal or boundary extraction. It provides the
evidence needed to evaluate future allele-selection rules without hiding
discarded internal SNP variants. Successful reruns with the option disabled
remove a stale table, and failed runs preserve the previous result bundle.

## Independent ±3 nt boundary candidates

Modern boundary inference now tests both junctions independently when a
database-assisted run or `--candidate-evidence true` forms its candidate set.
For each exact pre-collapse pair `R + S`, repeat start and repeat end each
move by −3…+3 bases: **49 combinations**, including the original cut.
The repeat becomes circular interval `[a, |R| + b)` and the spacer becomes
`[|R| + b, |R| + |S| + a)`. Each end can move three bases, so the total spacer
length can change by up to six bases. Equal shifts move the spacer without
changing its length. All bases come from the observed circle.

`--boundary-radius 0|1|2|3` controls this range (default 3).
`--boundary-policy catalogue` is the default additive selection policy when
`--array-database` and quality-bearing reads are available. Every candidate
is checked for exact repeat and spacer catalogue membership; promotion needs
the complete read context and same-source-array corroboration. The focal
spacer and every alternate cut of its circular core are excluded from its
corroborating anchors. The score is source category weight × confidence ×
`min(other distinct spacers, other distinct cores, 3) / 3`.

The selected policy requires a score of at least 0.5 (at least two other
known spacers/cores for Bona-fide evidence). It keeps an unambiguous highest
scoring cut in each nested ±1, ±2 and ±3 window, then takes their union.
Ties between different cuts within a window remain unresolved. Different
windows can prefer different cuts: these are labelled alternatives, not
verified boundaries. Groups require two distinct canonical spacers and obey
the explicit output length constraints. Existing calls and row-indexed views
are preserved; accepted additions appear in `boundary_supported` groups.

The other policies support reproducible ablations:

| Policy | Behavior |
|---|---|
| `off` | Disable promotion; candidate diagnostics still run when requested. |
| `catalogue-best` | One unambiguous catalogue winner over the entire requested range, with the weaker one-other-anchor threshold. |
| `catalogue-all` | Retain all eligible known same-array cuts, including competing boundaries. |
| `anchored-best` | One unambiguous winner; permit novel spacers anchored by other known same-array spacers. |
| `anchored` | All eligible known or anchored novel cuts. |

These alternatives are experimental. The broad read-supported lattice is
not emitted as calls: exact read presence alone cannot choose a biological
boundary. Missing reads/FASTQ qualities or a missing catalogue cannot promote
boundary candidates. Database inference with read validation disabled also
does not enable this selection stage automatically.

`boundary_evidence.tsv` links added rows to catalogue and read evidence.
With `--candidate-evidence true`, `boundary_candidates.tsv` records every
combination, its seed ID, shifts, sequences, discovery-length eligibility,
read support, independent repeat/spacer database scores and same-array scores.
Out-of-range combinations remain in this diagnostic table but cannot be
promoted (repeat 20–55 nt; spacer 20–70 nt). Seeds are sorted exact oriented
pre-collapse pairs and use the same numbering as `candidate_alleles.tsv`.
Different seeds can generate the same pair; evidence and selection deduplicate
exact pairs and circular cores rather than adding those duplicate votes.
Read counts cap at 32 and include saturation flags. Full contexts are shared
only where their literal sequences are identical; each valid variant receives
the same result as an independent exact-context scan.

Benchmark references and source genomes are used **only by the external
evaluation scripts**, never by the executable. The full development experiment
and reproducible policy comparisons are in
`ArraySearch/zuck_boundary_20261113` in the enclosing workspace.

On the frozen assembly-excluded development benchmark, the selected policy
raises exact CI recovery from 1,145 to 1,150 of 1,517 spacers and exact CCF
recovery from 1,592 to 1,604 of 1,703. It adds 20 distinct source-supported
spacer types, preserves all 6,605 previous rows, and improves exact reference
F1 in both panels. Strict recall also increases, while strict F1 falls
slightly. See `DEVELOPMENT.md` for the full tradeoff; these are development
results, not held-out biological accuracy estimates.

## Repeat, spacer and boundary confidence

Every run writes `element_confidence.tsv`, joined to `entries.tsv` by
`array_id` and `entry_index`. `--confidence-scores false` disables this
annotation and removes a stale confidence table on a successful rerun. All
original calls, IDs, and grouping views are retained. Cycle/graph replays use
`--grouping-reads` for read evidence; full runs use their input reads.
When a supported-call profile is enabled, its internal evidence is still
computed even if this annotation table is disabled. The two outputs share
one evidence pass. Use `--call-profile off --confidence-scores false` to
disable this pass; other enabled evidence features can still scan reads.

These values are **evidence indices, not calibrated probabilities**. They
answer different questions:

| Field | Interpretation |
|---|---|
| `repeat_score`, `spacer_score` | Quality-qualified exact sequence presence in reads; this alone does not establish a biological repeat or spacer. |
| `context_score` | Support for the entire literal repeat + spacer + first 16 bases of the next repeat, or its whole-span reverse complement. |
| `repeat_database_score`, `spacer_database_score` | Strongest eligible catalogue category-weighted source confidence for that exact canonical sequence; zero means no eligible evidence. |
| `repeat_family_score` | Saturating number of distinct read-supported circular cores sharing this exact repeat. The focal core is included; it is not an independent boundary vote. |
| `database_score` | Same-source-array corroboration from other distinct spacers and circular cores, excluding the focal spacer and core. |
| `boundary_score`, `boundary_status` | Catalogue support for the complete repeat/spacer decomposition, discounted when other emitted decompositions of the same circular sequence also have catalogue support. |

Read support uses distinct complete read sequences or ordered sequence pairs,
plus distinct start offsets. Identical sequence duplicates do not accumulate
votes. Each matched span requires all base qualities at least Phred 10 and
total expected errors at most one. The read index is
`min(sequence_units,4)/4 * min(start_offsets,2)/2`. Counts cap at 32; saturation
flags mean **at least** that many, not exact depth or independent molecules.
`NA` means the evidence was not assessable, including missing FASTQ qualities
or a span longer than every quality-bearing read. Zero means no qualifying
support under the stated rule, not proof of an incorrect call.

The catalogue index uses one complete source record, never accumulated votes
from duplicate assemblies. Its score is source category weight times source
confidence times `min(min(other_spacers,other_cores),3)/3`. Bona-fide evidence
requires at least one other spacer/core; Possible evidence requires at least
two. `boundary_score` divides this by the number of catalogue-supported
decompositions in the same circular class. Missing catalogue evidence leaves
unresolved boundaries unscored; a single emitted interpretation is not proof
that no competing cut exists.

Several shifted boundaries can have the same complete read context. Their
read evidence therefore ties even at high coverage. Boundary scores describe
the complete decomposition, not separate calibrated start/end probabilities.
`catalogue_preferred` means a catalogue preference, not a verified boundary.
The confidence experiment retains ambiguous calls because removing them lost
supported sequences and reference recoveries. Its scripts and full comparison
are in `ArraySearch/zuck_confidence_20261113` in the enclosing workspace.

## Supported calls and complete candidate evidence

`--call-profile standard` is the default additional output view. It selects
individual observed repeat/spacer rows whose spacers are 26–50 nt inclusive
and whose complete literal context has at least two quality-qualified
sequence units and two start offsets. The context and quality rules are the
same as the confidence annotation above. `--call-profile context` applies
the context rule without a profile length bound. `--call-profile off`
disables the four supported-call files and removes stale copies on a
successful rerun. The settings-file key is `call-profile`.

`--family-repeat-context true` is an experimental, disabled-by-default context
mode. The whole upstream repeat and spacer must still occur exactly in a
quality-qualified read. The following 16-base repeat anchor may additionally
come from an observed repeat allele in the same output family, after strand
orientation, with the same length and at most four substitutions. This handles
adjacent repeat copies that differ without accepting arbitrary mismatches or
borrowing alleles from other families. Read units and offsets are deduplicated
across alternative anchors. The mode changes final confidence and supported
selection only; discovery, inclusive rows, candidate evidence and database
validation retain their existing rules. Its policy is recorded in the run
manifest. Small controlled tests verify the mechanism; a broad benchmark
improvement has not been established.

All profiles preserve `CRISPR_Arrays_*.txt`, `arrays.tsv`, `entries.tsv`, and
the original grouping views. Explicit `--spacer-min-length` and other
output constraints still act upstream and can filter that inclusive output;
the call profile itself only adds a selected view. It does not expand the
candidate search, repair a boundary, or require a catalogue.

| File | Contents |
| --- | --- |
| `call_decisions.tsv` | Every emitted row, its original IDs, selection flag, context status and observability, bounded support counts, length eligibility, and decision reasons. |
| `supported_entries.tsv` | Selected oriented repeat/spacer rows, retaining original `array_id` and one-based `entry_index`. |
| `supported_arrays.tsv` | Selected membership within each original family, with consensus, entry count, distinct canonical spacer count, and explicit `singleton` or `multiple_spacers` status. |
| `supported_spacers.fasta` | Sorted unique reverse-complement-canonical spacer sequences supported by at least one selected row. |

The canonical FASTA does not validate every repeat pairing that contains the
same spacer. Family membership is a content view, not a verified genomic
locus or biological spacer order. Singleton families remain explicitly
labelled; no additional minimum family-size rule is applied. Enabled
experimental supplemental views can also have selected rows; join their
original IDs to `arrays.tsv` for `kind` and `parent_array_id` provenance.
Read-linked and database groups describe the inclusive rows and must be
intersected with selected original row IDs before evaluating selected-group
recovery. Do not transfer their original completeness scores to this view.

Missing reads, absent FASTQ qualities, or a full context longer than every
quality-bearing read produce `unknown`, with reasons `no_reads`,
`no_quality`, or `span_exceeds_max_quality_read`. Such rows remain in the
inclusive evidence and decision table but are not asserted as supported.
Observable spans lacking the required counts are `insufficient`; this does
not establish that a biological spacer is absent. Reasons also distinguish
`spacer_length_outside_profile`, `insufficient_context_units`, and
`insufficient_context_offsets`. A selected row has reason `passes_profile`.

The 26–50 nt profile is a transparent reporting domain, not a universal
biological spacer-length limit. The selected view can improve concordance
with a reference panel while excluding real, atypical, or sparsely covered
sequences. Preserve the candidate files for broader discovery and validate
the chosen profile on independent communities. Neither support indices nor
selection flags are calibrated correctness probabilities.

On the existing development community, a fresh 35,132,388-read native run
produces 1,757 supported spacer types and 1,586 exact union reference matches:
83.25% F1, versus 57.27% for the retained inclusive inventory and 75.49% for
the frozen CRASS comparison. The supported view loses three inclusive exact
matches and reduces strict clean-group recovery from 36/36 to 31/30 in the
CI/CCF panels. These are annotation-concordance results on reused development
data, not held-out biological accuracy. In a new mixed-length synthetic
fixture, the standard profile lowers F1 by excluding real 22- and 55-nt
spacers; `context` retains all recovered lengths. See `DEVELOPMENT.md` and
`ArraySearch/zuck_supported_calls_20261114/REPORT.md` in the enclosing workspace
for the full results, native outputs, and reproducible validation.

## CRISPRidentify array database

For source-linked annotations of **all** output repeats and spacers, and a
candidate viral protospacer search using CRISPRloci v2's CRISPRspacer engine,
see [Database matches and protospacers](docs/database-and-protospacer-matches.md).
The annotation stage works independently of the inference option below.

`--array-database arrays.tsv` enables confidence-weighted interpretation of
observed circular sequences. It adds labelled database-assisted groups and
retains all de novo groups, oriented rows, and read-link views. With this
option omitted, sequence calls are unchanged.

Prepare a sample-specific projection from an immutable CRISPRloci catalogue
snapshot, using cycles from a first Zuck run:

```bash
python3 scripts/build_array_database.py \
  --snapshot /path/to/crisprcasdb-20260927 \
  --cycles result/cycles.txt --output arrays.tsv --workers 4
build/Zuck --cycles result/cycles.txt --array-database arrays.tsv \
  --grouping-reads reads_R1.fastq reads_R2.fastq \
  --output-folder database-result --threads 8 --ram 4G
```

The Python preparer needs Python 3.11 or later with SQLite deserialization
support. It streams and verifies every compressed source shard, selects arrays
using both strands of sample-observed circular 20-mers, and writes a bounded
TSV plus a SHA-256 provenance manifest. No benchmark annotation is used to
select the projection. The native inference checks complete sequences; a
20-mer seed is insufficient for a call. Preparation is a separate,
sample-dependent cost. A prepared TSV also works with `--input-files` and
`--graph`, but a projection prepared for a different sample can omit relevant
arrays.

An accepted boundary must explain a complete observed cycle using an exact
catalogue repeat and spacer that occur in the same source array. By default,
only Bona-fide arrays with confidence at least 0.75 and at least two distinct
matched spacers can guide a call. `--database-include-possible true` also admits
Possible arrays with confidence at least 0.5 and at least three distinct
matches. These source confidence scores are heuristic, not probabilities.

Each candidate first needs exact FASTQ support for the full repeat, full
spacer, and the first 16 bases of the next repeat, or the reverse complement
of that entire span. All span bases must have Phred quality at least 10 and
their expected errors must sum to at most one. Support requires at least two
distinct full read or read-pair sequences and two start offsets. These units
are not proven independent molecules. FASTA cannot supply this quality
evidence. Cycle replays therefore need `--grouping-reads`; full read runs use
their input reads automatically. `--database-read-validation false` is an
explicit diagnostic mode and does not emit database content views.

After read validation, repeat alleles and known spacers contribute separately
to a bounded evidence score:
`confidence × category_weight × (2 × min(spacer_matches, 4)/4 + min(repeat_matches, 2)/2)`.
Category weight is 1 for Bona-fide and 0.5 for Possible. Counts refer to distinct
canonical sequences within one source array. Duplicate assemblies and repeated
cycle traversals cannot multiply the score. Equal-scoring competing boundaries
are left unresolved. A catalogue spacer absent from the observed cycles is
never copied into output.

`--database-novel-spacers true` is an experimental arm for observed spacers
absent from a supporting source array. It requires the same exact read context
and a repeat with the required known same-array spacer anchors. A novel spacer
cannot count as its own known anchor; eligible known interpretations take
precedence. It is disabled by default.

`database_evidence.tsv` identifies the source-array signature, confidence,
distinct repeat/spacer match counts, score, cycle support, and read support for
each assisted row.
The database input and effective thresholds are recorded in the run manifest.
These groups provide sequence content evidence, not verified genomic loci or
spacer order. The supplied September snapshot contains accepted CRISPRidentify
2.3.2 region-screen arrays; it is distinct from the later whole-genome campaign.

Database content views use same-array membership within a read-validated
exact-repeat parent. They require at least three matched spacers for Bona-fide
sources or four for explicitly enabled Possible sources. Equal sets are
deduplicated, nested subsets and full-parent duplicates are omitted, and a
parent with more than 16 incomparable maximal sets is left unresolved. The
`database_groups.tsv` and `database_group_members.tsv` sidecars reference
existing rows; they introduce no new sequences. Disable these views with
`--database-groups false`.

For an assembly-excluded evaluation, provide `--exclusions exclusions.json`
and `--excluded-output arrays.excluded.tsv` to the preparer. The JSON must
contain `exclusion_accession_roots`, including both RefSeq and paired GenBank
accession roots. Exclusion occurs before deduplication. Related genomes can
still contribute shared repeat/spacer sequences.

## Results and reproducibility

Every successful run, including a zero-array run, writes:

- `CRISPR_Arrays_N.txt`: grouped repeat/spacer rows and consensus sequences.
- `arrays.tsv`: one row per output group, with counts, consensus, kind, and
  parent identifiers for experimental alternative views.
- `entries.tsv`: the exact repeat variant and spacer for every output row.
- `run_manifest.json`: source fingerprint, compiler, effective settings,
  stage/filter counts, input-cycle checksum, and output inventory/checksums.

When read-supported grouping runs, `read_linked_groups.tsv` references the
original parent rows and `read_linked_edges.tsv` records the supporting
sequence-unit counts, junction-start diversity, and quality rule. These
sidecars do not add spacer discoveries or rewrite the primary tables.

Table identifiers match the numbered text files. Tables are deterministic
for a fixed replay input and configuration; timestamps and timing vary.
The compiled source identity is SHA-256 over sorted project/vendored source
paths and their digests. It includes local changes and refreshes on rebuild,
but excludes externally fetched dependency trees. File checksums use XXH64
with seed 0 for accidental-change detection.

Replay manifests identify the cycles actually processed and mark graph-search
settings as not applied in that run. Read and external graph paths are
recorded without hashing those large inputs. The built-in `scoring` field
describes confidence annotation (or is null when disabled); record external
benchmark scores separately. Stage observations can
include redundant graph traversals; they are not read depth.

Results are staged before replacing an existing bundle. Parsing and staged
write failures preserve earlier results, and rename failures trigger rollback.
Successful reruns remove stale numbered files, preserve unrelated files, and
reject symlink/directory destinations. Use one writer per output directory:
replacement across multiple files is not atomic for concurrent readers or
process/power failure.

## Interpretation and limits

Groups describe sequence content. They do not establish unique genomic loci,
terminal repeats, biological spacer order, or calibrated confidence. Output
repeat rows are cycle-derived evidence, not observed genomic repeat copies.
When redundant strand-equivalent group boundaries are coalesced, all observed
oriented repeat/spacer rows are retained. A resulting group can contain both
observed orientations; its rows remain unordered sequence evidence.
Tandem-repeat backgrounds and small-family errors remain possible. Search depends on logical partition count and visited history; worker count
and effective OpenMP team size do not change those histories.

The 12 November benchmark used one development community and incomplete
computational reference panels. Its concordance/F1 scores are not biological
precision or evidence of performance on unseen communities. See
`DEVELOPMENT.md` for the measured baseline and research priorities, and
`UPSTREAM.md` for authorship and dependency acknowledgements.

```bash
docker build -t zuck .
docker run --rm --user "$(id -u):$(id -g)" \
  -v "$PWD/data:/data:ro" -v "$PWD/result:/out" \
  zuck --input-files /data/reads_R1.fastq /data/reads_R2.fastq \
  --output-folder /out --threads 8 --ram 4G
```

The final image is `scratch`-based and contains the same static executable as
the standalone export. It intentionally has no shell or package manager.
