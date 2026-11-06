# my_mcaat development log

This is the single living record for council decisions, implementation
iterations, and benchmark outcomes. Disposable builds and runs stay outside
the repository under `/tmp`.

**Correction to the previous vNext interpretation:** the experimental
alternative-repeat sidecar did not discover four new spacers. They and the
observed repeat allele were already present in the primary output; the old
comparison code had skipped that allele because its consensus failed a
preliminary similarity gate. With allele-aware strict scoring, the primary
call is split/merged and the sidecar improves the same reference's grouping to
clean without changing its supported-spacer count. The historical scores
below document the old comparison, not an independent sensitivity gain.

## Frozen baseline

Input: the existing pooled InSilicoSeq paired-read metagenome (45 assemblies,
99 replicons, 17,566,194 read pairs, 20–40× assigned coverage).

MCAAT v1.0.0 produced 413 blocks and 5,098 spacer rows in 1,406.2 seconds.
The frozen automated scorecard reports:

| Stratum | Clean | Split/merged | Partial | Miss | Supported spacers |
|---|---:|---:|---:|---:|---:|
| CRISPRidentify Bona-fide | 15 | 16 | 21 | 6 | 1,261 / 1,560 |
| CRISPRCasFinder Evidence-4 | 18 | 12 | 27 | 2 | 1,429 / 1,741 |

Council adjudication remains separate: CRISPRidentify strict has 4 faithful
and boundary rescue has 7; CRISPRCasFinder Evidence-4 has 22 faithful, 11
merged, 1 split, and 25 under-90 partial arrays. The automated and adjudicated
numbers must not be mixed.

## Council round 1 — 2026-11-05

### Finding

The default 77-nt maximum cycle length is a hard sensitivity ceiling. All
207/207 CRISPRCasFinder Evidence-4 spacer elements with repeat+spacer length
above 77 are complete misses; 221/222 such CRISPRidentify Bona-fide elements
are misses. A limit of 90 makes 206/207 and 208/222 of those elements,
respectively, geometrically reachable. A limit of 100 adds little expected
recall for substantially more branching risk.

The cycle enumerator also erased all accumulated candidates when a start node
reached its 500-cycle safety cap. Saturation must retain the capped evidence.

### Iteration 1 plan

1. Create the independent `my_mcaat` source snapshot and preserve exact
   upstream revisions.
2. Add a cycles-only replay mode; retain cycles and remove graphs by default.
3. Raise the default maximum cycle length from 77 to 90.
4. Stop at 500 cycles per start node without clearing retained cycles.
5. Parameterize the frozen comparison program and write only a compact
   scorecard during iterations.
6. Run the full pooled benchmark once, then replay later post-processing
   experiments from its retained cycles.

### Iteration 1 acceptance gate

- at least 1,570/1,741 council-supported Evidence-4 spacers (the automated
  scorecard baseline is 1,429 and is reported separately);
- at least four current Evidence-4 partial arrays cross 90% union recall;
- at least 25/59 council-faithful Evidence-4 arrays;
- no increase above 11 merged and 1 split council calls;
- runtime and peak memory no more than twice baseline;
- identical settings and thread count produce stable array/spacer tables.

### Deferred risks

- global shared-23-mer spacer union can destroy locus-specific evidence;
- repeat edit-distance single linkage can chain distinct families;
- block-wide suffix trimming can harm already-correct boundaries;
- per-thread visited state makes candidate suppression order-dependent;
- faithful precision and biological order require synthetic truth and later
  streaming read/path provenance, not only caller agreement.

## Iteration history

| Iteration | Change | Runtime | Automated CI BF | Automated CCF E4 | Decision |
|---|---|---:|---|---|---|
| 0 | MCAAT v1.0.0 baseline | 1,406.2 s | 15 clean; 1,261/1,560 spacers | 18 clean; 1,429/1,741 spacers | frozen |
| 1 | max cycle 90; retain capped cycles; replay support | 1,424.6 s | 19 clean; 1,360/1,560 spacers | 22 clean; 1,522/1,741 spacers | partial gain; gate not met |
| 2d | exact support aggregation + family-local 23-mer collapse | 6.2 s replay | 19 clean; 1,366/1,560 spacers | 22 clean; 1,531/1,741 spacers | modest gain; gate not met |
| 3a | support-ordered complete-link repeat families | 6.2 s replay | 22 clean; 1,439/1,560 spacers | 24 clean; 1,603/1,741 spacers | substantial gain; one clean call short of gate |
| 3c | support-anchored exact-containment boundary correction | 6.0 s replay | 24 clean; 1,472/1,560 spacers | 25 clean; 1,635/1,741 spacers | gate passed; boundary safety audit required |
| 3e | terminal-only correction, row-specific trimming, support preservation, consensus coalescing | 1,756.7 s full; 6.1–13.1 s replay | 23 clean; 1,438/1,560 spacers | 25 clean; 1,603/1,741 spacers | accepted |

## Iteration 1 outcome

The graph was identical to baseline (634,337,564 nodes, k=23; 297,839,652
valid edges after filtering). `my_mcaat` recovered 233,877 raw cycles across
8,078 start nodes, versus 112,796 cycles across 7,318 starts in MCAAT. Runtime
rose by only 18.4 seconds (1.3%). The final output contained 429 arrays and
4,957 spacers.

Automated high-confidence results improved materially but missed the iteration
gate:

- CRISPRidentify Bona-fide clean arrays increased 15→19 and supported spacers
  1,261→1,360 (+99; 80.8%→87.2%).
- CRISPRCasFinder Evidence-4 clean arrays increased 18→22 and supported spacers
  1,429→1,522 (+93; 82.1%→87.4%).
- Evidence-4 split/merged calls increased 12→13 and complete misses increased
  2→3, so structural errors did not remain flat.

Despite 2.07× more raw cycles, final spacer rows fell from 5,098 to 4,957.
This strongly supports the predicted destructive interaction between new
candidate volume and the global shared-23-mer union-find. Iteration 2 should
use cycles-only replay to isolate post-processing changes without rebuilding
the graph.

## Council round 2 and iteration 2 — 2026-11-05

### Finding

The global spacer union assigned each transitive 23-mer component to only one
repeat variant. That could discard the same spacer evidence from compatible or
unrelated repeat families. It also had a measured blind spot: all 708 exact
within-array duplicate rows in iteration 1 were shorter than 23 nt and could
not enter the deduplication index.

The council therefore requested exact support aggregation by repeat variant,
family-scoped spacer handling, deterministic tie-breaking, and a second exact
collapse after repeat-tail trimming. It also recommended testing bounded
repeat clustering separately from spacer handling.

### Replay experiments

All disposable outputs were kept under `/tmp`; only this summary is retained.

1. Exact-only family aggregation exposed 125,388 exact repeat/spacer alleles
   from 225,783 observations and emitted 118,664 rows. This was rejected as an
   unusable upper-recall/noise bound before promotion.
2. Conservative complete-link near-duplicate clustering (>=90% reciprocal
   span, >=95% overlap identity, <=2 mismatches) still emitted 51,328 rows and
   was rejected as insufficiently selective.
3. Scoping the established shared-23-mer collapse to each compatible repeat
   family emitted 556 arrays and 5,331 spacers in 6.2 seconds. It preserved
   family evidence and removed exact short duplicates, but improved the frozen
   scorecard only to 1,366 CI Bona-fide and 1,531 CCF Evidence-4 supported
   spacers. Clean counts stayed 19 and 22; CCF split/merged rose 13→14 and CI
   complete misses rose 5→6. This is a correctness improvement, not yet a
   substantial detection improvement.

The comparison program now accepts either one MCAAT-compatible output file or
a directory of numbered output files, allowing large stress-test candidates
to be parsed without creating concatenated artifacts.

## Council round 3 and iteration 3 — 2026-11-05

### Finding

Single-link repeat clustering was the next destructive step: 45/556 repeat
families had a member-to-member edit-distance diameter above four, with a
maximum of fourteen. Retained cycles also contained exact evidence for
130/210 still-missing Evidence-4 spacers. Several of those rows had a repeat
prefix or suffix assigned to the spacer because the cycle was cut in a
different repeat phase.

The council replaced repeat single linkage with deterministic, support-ordered
complete linkage. A repeat joins a family only when it is within edit distance
four of every existing member. Ties are resolved by maximum distance,
aggregate support, and lexical anchor. This made the diameter bound an
invariant and raised iteration 3a to 1,439 CI Bona-fide and 1,603 CCF
Evidence-4 supported spacers.

Iteration 3b tested exact-containment-only spacer collapse and emitted 1,045
arrays with 67,568 spacer rows. It was rejected as an unselective upper bound.
Iteration 3c restored family-local shared-23-mer collapse and used a
best-supported allele plus an at-most-eight-base containment correction. It
crossed the provisional score gate, but council audit found that unrestricted
substring matching and block-wide trimming could over-credit internal matches
or alter rows that did not support the inferred tail.

### Frozen substantial-improvement gate

- CCF Evidence-4: at least 1,570 supported spacers and 25 clean arrays.
- CI Bona-fide: at least 1,404 supported spacers and 22 clean arrays.
- Every raw repeat family has edit-distance diameter at most four.
- No within-array exact duplicate spacers; deterministic normalized output.
- No post-boundary-correction spacer shorter than the 20-nt extraction floor.

## Council round 4 and accepted iteration 3e — 2026-11-06

The final safety pass made every shortening explicit and local:

1. Equal-support selection starts from the longest observation; a shorter
   allele may replace it only when it is an exact prefix or suffix (or the
   reverse-complement equivalent) and differs by at most eight bases.
2. A validated family tail is removed only from rows whose own suffix matches
   that tail within edit distance two, and original observation support is
   preserved through trimming.
3. Boundary correction cannot create a spacer shorter than 20 nt.
4. Raw families that converge to the same final repeat consensus are
   coalesced deterministically with exact spacer deduplication.
5. Reused output directories remove only numbered `CRISPR_Arrays_N.txt`
   files, including on a valid zero-detection replay.

The accepted retained-cycle replay contains 661 arrays and 7,212 spacers. Its
frozen score is:

| Stratum | Clean | Split/merged | Partial | Miss | Supported spacers |
|---|---:|---:|---:|---:|---:|
| CRISPRidentify Bona-fide | 23 | 16 | 13 | 6 | 1,438 / 1,560 |
| CRISPRCasFinder Evidence-4 | 25 | 13 | 18 | 3 | 1,603 / 1,741 |

Relative to upstream MCAAT this is +8 and +7 clean arrays, and +177 and +174
supported spacers, respectively. The final validation deliberately discarded
the extra low-span exact-containment credit seen in 3c: the accepted gain is
therefore smaller but has stronger boundary semantics. The substantial-
improvement gate still passes in both strata.

Two independent final replays were identical after removing only the generated
timestamp header. Reproducibility identifiers:

- retained cycles SHA-256:
  `b1bcb5cf08fcdda3f41b511aeb169c68a8a18304049d6093c1ac1da897ce2a57`;
- normalized final output SHA-256:
  `8c905ac91bb301c813003d2b3eb92959eb2d56fc27653560067e0b8a102be6c6`;
- final release image/binary build:
  `sha256:013c9b19e53543e2e3edaeb64cc5a40916b7aeee760b32ae0e5f5e41000199d7`.

Structural checks found no malformed blocks or within-array exact duplicates;
emitted spacer lengths are 20–67 nt. A zero-cycle hygiene test removed stale
numbered results while leaving an unrelated similarly named file untouched.

### Independent end-to-end robustness run

A second run started from the original 17,566,194 paired reads with 10 threads,
an 8 GB requested RAM budget, a 12 GB container limit, and maximum cycle length
90. The graph exactly reproduced the baseline topology (634,337,564 nodes at
k=23; 297,839,652 valid edges after filtering). Total runtime was 1,756.7
seconds, 1.25× the upstream MCAAT baseline and below the 2× ceiling. The run
overlapped compilation and scoring jobs, so this is a conservative wall-clock
measurement rather than an isolated speed benchmark.

Parallel cycle enumeration is not byte-deterministic: this run found 240,034
cycles across 8,067 start nodes instead of the retained run's 233,877 cycles
across 8,078 starts. This confirms the deferred per-thread visited-state risk.
The accepted post-processing is deterministic for a fixed cycle file and the
second cycle set independently passed the improvement gate:

| Stratum | Clean | Split/merged | Partial | Miss | Supported spacers |
|---|---:|---:|---:|---:|---:|
| CRISPRidentify Bona-fide | 24 | 16 | 12 | 6 | 1,473 / 1,560 |
| CRISPRCasFinder Evidence-4 | 26 | 13 | 17 | 3 | 1,637 / 1,741 |

It emitted 645 arrays and 7,660 spacers, again with a 20-nt minimum and no
within-array exact duplicates. The independent cycle and normalized output
SHA-256 values are, respectively,
`878a4c52625d9dc4fb55ca178e8c091740da7a098ca36f88d4c72739a4ec469e`
and `8eb1bb9575196351e68dcda41921d73aa1a9515511070d939600adeffae22eed`.
Because both full runs exceed the gate but raw cycle sets differ, the retained
first run remains the frozen promoted artifact; the second is a robustness
check, not a best-of-two replacement.

## vNext council round 1 — 2026-11-06

The v1.0.0 release commit `0650fe5` remains frozen. Development continues on
`performance-next`, and disposable builds, graphs, scorecards, and candidates
remain under `/tmp`.

Four independent council tracks reviewed detection algorithms, systems,
release-error diagnostics, and validation. Their combined findings were:

1. **Graph memory was passed in the wrong unit.** `Settings::ram` is GiB, but
   MEGAHIT's `host_mem` is bytes. An `--ram 8G` run therefore passed roughly
   eight bytes and silently fell back to a minimum-memory sorting plan. Graph
   build/load consumes about 94% of cold runtime, so checked GiB-to-bytes
   conversion plus cgroup-aware defaults is the highest-value systems fix.
2. **Cycle storage is 94.6% redundant.** All 17,978,325 adjacent k-mer pairs
   in the retained release evidence overlap by 22/23 bases. One glued sequence
   per cycle is lossless and is estimated to reduce 437,092,848 bytes to about
   23.6 MB. The new writer must remain readable alongside legacy lines and
   reproduce the same normalized array output.
3. **Traversal duplicates are not independent support.** The 233,877 release
   lines reduce to 46,259 unique circular sequences after rotation and
   reverse-complement normalization. Canonical evidence/provenance is a later
   foundation, but naïvely coalescing reverse-complement output families was
   rejected because it reduced CI clean calls from 23 to 22.
4. **The current physical-thread visited bitmap is both large and semantic.**
   It occupies about 756 MiB at ten threads on the frozen graph and makes
   retained evidence schedule-dependent. Deterministic ownership, bounded-BFS
   prefiltering, and diversity-aware top-K retention are promising graph-stage
   experiments, but they must remain isolated detection changes.
5. **Blind shared-23-mer transitivity is the safest first replay target.** A
   shared k-mer should propose a spacer pair, not prove that every member of
   an A–B–C chain is one allele. The selected experiment uses support-ordered
   star or complete-link verification with full-span edit/coverage checks and
   explicit terminal-containment handling. The exact v1.0.0 transitive mode is
   retained behind an ablation switch.
6. **Boundary inference has high upside and high gaming risk.** Robust
   two-sided repeat/spacer change-point inference could rescue phase-shifted
   evidence, but must be evaluated after clustering and only with strict
   reciprocal-length matches, row-local provenance, and synthetic boundary
   truth. It is not bundled into iteration 1.

The release-error trace found 260 missing high-confidence spacer occurrences
(122 CI Bona-fide and 138 CCF Evidence-4). Of these, 151 had no exact retained
cycle occurrence, 11 were lost during rotation/grouping/extraction, and 98
survived extraction but disappeared downstream. `NC_010175_3` alone has 38
post-extraction losses and is the primary replay diagnostic. In contrast,
`NC_009495_4` has only one of five exact raw-cycle spacers and
`NZ_CP076401_7` has none, so replay-only tuning must not claim to solve those
raw-evidence ceilings.

### Iteration 1 plan and gates

- fix the RAM unit and cgroup default in isolation;
- write deterministic compact cycle records and count cap-saturated starts;
- compare `transitive`, support-star, and verified complete-link spacer
  components from the identical frozen cycle file;
- keep all release floors: CI 23 clean / at most 16 split-merged / at most 6
  misses / at least 1,438 supported; CCF 25 / 13 / 3 / 1,603;
- reject any clean-to-partial/miss regression, spacer below 20 nt, duplicate
  row, repeat-family diameter violation, or uncontrolled output growth;
- require either a new clean high-confidence array or at least one percentage
  point strict best-array F1 improvement without regression in the other
  stratum; union recall or exact-containment-only gains do not qualify;
- for speed, require identical graph topology and detection output, with at
  least 15% isolated median end-to-end improvement or 20% stage improvement.

Deferred techniques remain recorded for later council rounds: robust
multi-phase/two-sided boundary inference, cyclic repeat clustering, canonical
cycle support, deterministic independent enumeration, bounded BFS, adaptive
90–110 continuation, repeat-anchored bubble detection, read-pair junction
provenance, locus separation, and calibrated candidate ranking on held-out
truth rather than caller agreement.

## vNext iterations 1–2 and council round 3 — 2026-11-06

Iteration 1 validated the systems changes and rejected the first detection
idea. A fresh release-image graph run on the 17,566,194 paired reads produced
the expected 634,337,564-node graph and 297,839,652 valid post-pruning edges.
Its graph build/load took 1,363.2 seconds and its cycle stage took 99.4
seconds. Replaying that retained graph with the vNext writer took 97.4 seconds
for cycle search and wrote 238,316 records in a 23 MB compact file instead of
the release writer's 422 MB file. Replaying that compact file reproduced the
same biological output byte-for-byte after excluding the generated timestamp.

The compact writer is headerless for v1.0.0 parser compatibility and verifies
the 22-base overlap before using one-base continuation encoding. Any future
violating path falls back to its lossless legacy k-mer line. The automatic RAM
default now uses `MemAvailable` and the process' nested cgroup remaining
budget; explicit GiB values are converted to bytes before entering MEGAHIT.

Complete-link and support-star component splitting were rejected before
scoring because they emitted 51,583 and 42,165 spacers respectively, compared
with the release's 7,212. They violated the output-growth gate by several-fold
and took 32–43 seconds per replay. The default was therefore restored to exact
release-compatible transitive components.

Council round 2 proposed a bounded additive experiment: preserve all release
primaries, repeats, tail decisions, and accepted arrays, then permit only
support-backed satellite medoids under global, family, component, and array
caps. The final iteration-2 candidate was structurally safe:

- 661 arrays and 7,428 spacers;
- all 7,212 release rows retained;
- 216 additions, exactly the 3% ceiling;
- no novel repeat fields, duplicates, sub-20-nt spacers, or local cap breach;
- 6.9-second frozen-cycle replay.

The biological gate nevertheless rejected it. CI Bona-fide support moved only
from 1,438 to 1,440 and CCF Evidence-4 from 1,603 to 1,605, with no new clean
array. Macro best-array F1 rose only about 0.00117 in each stratum, CCF purity
q10 declined, and only three new canonical high-confidence spacer sequences
were associated with 216 additions (1.39% marginal precision). With the
exact-containment override disabled, the candidate gained only two CI and one
CCF supported spacers and no clean status.

The apparent `NC_010175_3` best-array gain of eight was a boundary artifact:
all eight new rows ended in `TTTCA`, which matches the repeat in the opposite
orientation within one edit. They reached only 0.872–0.885 reciprocal length
coverage and none survived strict reciprocal-span matching. More broadly,
132/216 additions retained a 5–16 nt repeat-like terminal suffix. This evidence
rules out selecting a smaller budget after seeing reference labels; that would
overfit the benchmark rather than improve detection.

Council round 3 therefore pivoted from extra representative emission to the
upstream boundary representation:

1. Rotate the true circular core, not the already closed linear record, and
   rebuild the 22-base closure after rotation.
2. Score that top-1 normalization in isolation before adding boundary logic.
3. If needed, infer repeat boundaries on both sides of the primary anchor with
   distinct canonical-cycle support and orientation-aware validation.
4. Only then consider a targeted second anchor; deduplicate hypotheses by
   source cycle so phases cannot inflate support.
5. Treat exact-containment-only gains as non-qualifying and require reciprocal
   identity/coverage, stable purity, and independent-cycle confirmation.

The medoid lane was rejected and its dormant implementation was removed before
promotion. Any future redesign requires canonical cycle provenance and a new
held-out evaluation; none of that rejected lane remains in the active code.

## vNext iterations 3–4 and council round 5 — 2026-11-06

Iteration 3 tested the council's corrected circular-core normalization as a
global replacement. It was rejected. The output contracted from 661 arrays /
7,212 spacers to 513 / 4,761. Containment-enabled scoring appeared to gain
some supported spacers, but strict reciprocal-span scoring collapsed from
15 to zero clean CI Bona-fide arrays and from 17 to three clean CCF
Evidence-4 arrays. This established that corrected rotation contains useful
evidence but cannot safely replace the frozen primary boundary convention.

The failure audit exposed one narrow, label-free alternative view. A later
scorer audit established that the parent already contained this evidence; see
the correction at the top of this log. Release `Array_274`
contains nine accepted spacers and an exact observed non-consensus repeat
allele associated with four of them:

```text
parent consensus: CGTTTCTAATGTACCGTAGTGGAGTTGAAAC
observed allele:  GTTTCTAATGTACCGTGTGGAGTTAAAAC
```

The allele has scorer-equivalent oriented cyclic identity `0.689655` to the
parent consensus and reciprocal length coverage `0.935484`. Because identity
is below the existing `0.80` repeat-candidate gate, the parent consensus
cannot expose this already accepted four-spacer subarray to downstream
matching.

Iteration 4 therefore keeps the complete release result frozen and, only when
`--experimental-alternative-views` is supplied, emits a supplemental sidecar
when all of the following reference-independent rules hold:

- the candidate is an exact observed non-consensus repeat allele owned by one
  accepted parent after canonical orientation checks;
- it has at least four exact/RC-unique accepted parent spacers and is a strict
  subset covering between one third and two thirds of the parent;
- cyclic identity to the parent consensus is below `0.80`, reciprocal length
  coverage is at least `0.90`, and every spacer is at least 20 nt;
- no novel repeat or spacer sequence is created;
- global growth is capped at 1% of primary arrays and 0.25% of primary rows.

Supplemental views are written to the next numbered result file, never into a
primary file, and carry an experimental `Kind` header. Normal runs leave this
lane disabled, and the primary-only `getSystems()` interface excludes it. Two
opt-in replays at one and ten threads preserved the release primary
normalized SHA-256 exactly:
`8c905ac91bb301c813003d2b3eb92959eb2d56fc27653560067e0b8a102be6c6`.
The final supplemental file has normalized SHA-256
`204434d50b54f38d797bb4b12671ed97dbf60c849e143c26fd60be60c06312d2`.
It contains one `Array_662` with four rows, so the union is 662 arrays and
7,216 spacers. Frozen-cycle replay took 6.2–6.3 seconds. The verified image is
`sha256:3286b38dcb6ca01968481653db3f813f1410ee4982badf082615ad228b4c74c2`.
The supplemental digest differs from the pre-promotion candidate only because
the final sidecar now carries its explicit experimental `Kind` header. Two
final opt-in replays at one and ten threads produced identical normalized
primary and sidecar digests; a separate default-off replay emitted only the
release-identical primary file.

The standard high-confidence score changes only for the CCF Evidence-4 view
of the recovered locus:

| Stratum | Release | Iteration 4 | Supported spacers |
|---|---:|---:|---:|
| CI Bona-fide | 23 / 16 / 13 / 6 | 23 / 16 / 13 / 6 | 1,438 → 1,438 |
| CCF Evidence-4 | 25 / 13 / 18 / 3 | 26 / 12 / 18 / 3 | 1,603 → 1,603 |

Under the historical consensus-gated strict scorer, the CCF Evidence-4 vector changes from
`17 / 7 / 15 / 20` to `18 / 7 / 15 / 19`, with supported spacers increasing
from 1,043 to 1,047. CRISPRidentify independently annotates the same locus as
a Possible candidate and also changes it from miss to clean under strict
scoring. These are two caller confirmations of one biological locus, not two
independent discoveries. All four matched spacers have identity and
reciprocal coverage `1.0`; no reference status regresses.
The allele-aware audit later showed that the parent already contained those
four sequences, so the historical supported-spacer increase is a scorer
artifact, not newly detected sequence evidence.

The promotion claim is deliberately limited: this implementation presents
one previously hidden repeat-allele subarray without changing primary
detection or adding spacer evidence. The thresholds were frozen after the in-sample audit, but generic
sensitivity improvement still requires held-out genomes with positive and
negative allele-partition truth. The unrelated `NC_003210_1` best-repeat
identity changes from `0.6154` to `0.6552` without crossing any evidence or
status threshold and is recorded as a secondary association, not a gain.
The final council therefore approved the lane only as a default-off
experimental sidecar and stopped further tuning on this benchmark to avoid
converting one in-sample recovery into an overfit production default.

## Literature and error-audit council — 2026-11-06

### What changed in the interpretation

The comparison program considered observed repeat variants only after the
primary consensus passed a preliminary similarity gate. This concealed the
existing observed repeat and four existing spacer rows in `Array_274`. Under
strict allele-aware scoring the primary already supports those four spacers;
the experimental sidecar changes their reference from split/merged to clean
by presenting a smaller subset, with no new spacer evidence. The earlier
strict gain in supported spacers was a scoring artifact.

An exhaustive audit of 117 high-confidence references against all 662
candidate arrays found four reference/array pairs where the consensus gate
hid a repeat-supporting allele. Only two pairs had at least two spacer matches
under the permissive matcher; only one of those survived strict matching. The
15-mer candidate index can theoretically omit valid approximate matches, but
an exhaustive scan found no such omitted spacer matches among those 77,454
reference/array pairs. The allele-aware matcher must rank a qualifying
full-span repeat above a higher-identity truncated one; two pairs in this
corpus exposed that distinction, although neither affected array statuses.

The comparison workspace now has opt-in `--allele-aware` and `--strict-span`
flags; default behavior preserves historical scorecards. Strict allele-aware
primary results are CI Bona-fide 15/9/14/20 with 957 supported spacers, and
CCF Evidence-4 17/8/15/19 with 1,047. These vectors are clean/split/partial/
miss. The corresponding sidecar vectors are 15/9/14/20 and 18/7/15/19;
supported counts stay 957 and 1,047. The high-confidence reference callers
overlap substantially: 44/58 CI rows overlap 45/59 CCF rows on their source
contigs, so their results are not independent biological replicates.

### Literature and hypotheses

- [metaCRISPR](https://academic.oup.com/bioinformatics/article/32/17/i520/2450786)
  recruits reads and mates around repeats, clusters using long overlaps and
  pairs, and performs local overlap-graph assembly. It identifies shared
  repeats and heterogeneous abundance as significant ambiguity sources.
- [Crass](https://doi.org/10.1093/nar/gkt183) retains spacer graphs and
  explains how sequencing errors, truncated reads and shared sequence can
  create false branches; co-occurrence alone is not proof of adjacency.
- [MetaCRAST](https://pmc.ncbi.nlm.nih.gov/articles/PMC5592083/) mines reads
  given repeat seeds and is therefore complementary, but cannot detect an
  unseeded family.
- [MCAAT](https://pmc.ncbi.nlm.nih.gov/articles/PMC12342471/) uses repeat-
  anchored graph cycles and multiplicity, which combines copy number and
  coverage; thresholds can lose low-abundance or divergent evidence.
- [CRISPRidentify](https://pmc.ncbi.nlm.nih.gov/articles/PMC7913763/) reports
  false-positive repeat structures and demonstrates why precision controls
  cannot be replaced by sensitivity against another caller.

### Error trace and read-based experiments

The original error trace counted 260 high-confidence caller spacer
occurrences missed by the standard scorer: 151 lacked an exact retained
cycle, 11 were lost during rotation/grouping/extraction, and 98 survived
extraction but disappeared downstream. Of these, 246 were at least 20 nt and
eligible for a reference-guided read scan. All 246 had exact full-spacer read
matches and 243 also had a read containing both expected 8-base flanks.
These matches establish input sequence presence, not unique locus assignment;
the two callers overlap and the other 14 short CI occurrences were excluded.
All 1,618 CI repeats and 1,560 CI spacers, and all 1,800 CCF repeats and
1,741 CCF spacers, agreed with their source FASTA coordinates. That verifies
extraction, not the biological correctness of each caller's boundaries.

A diagnostic read-recruitment prototype used only uniquely owned primary
consensus repeats as seeds. It required exact repeat-spacer-repeat in a read,
quality >=20, spacer length 20–67 nt, at least two even-indexed read pairs and
one odd-indexed pair; a pair was counted once. Its 49
nonidentical parent-spacer candidates from 28 parents all occur as exact
repeat-spacer-repeat sequence in the simulated source genomes. Seven contain
a whole additional repeat inside the proposed spacer, 15 approximately
match an existing parent spacer under 90% identity/90% span, and three match
a high-confidence reference spacer with compatible repeat. The full candidate
changed neither caller's corrected strict status vector or supported count.
The split checks observation consistency, not independent validation: primary
seeds were already derived from the complete read set and the odd half also
participated in candidate selection. Exact genomic occurrence does not
establish a CRISPR boundary.
This prototype was not promoted.

A second pre-specified variant included uniquely owned observed repeat alleles
as seeds, retaining the same exact-read, base-quality and even/odd pair
support rules. It produced 184 canonical-distinct new parent-spacer rows in
60 parents, including 137 not present in the consensus-only candidate and
47 shared with it; two consensus-only rows were excluded by the stricter
cross-parent allele ownership rule. All 184 exact repeat-spacer-repeat
strings occur in the simulated source genomes, but 20 proposed spacers
contain an entire repeat, and 97 approximately match an existing parent
spacer at 90% identity and reciprocal span. Three match a high-confidence
reference spacer exactly and 15 approximately match one with a compatible
reference repeat. Neither high-confidence caller's strict allele-aware
status vector or supported-spacer total changed. The three exact matches may
already be supported elsewhere; these counts are not a precision estimate.
All 184 are exact/RC-distinct from the primary spacer output, but 83 only
meet the minimum two-plus-one read-pair support. The scratch parser assumes
sequential primary array IDs, which holds here but is not a general interface
guarantee. This variant was also rejected for promotion.

A separate same-read co-occurrence experiment used exact existing spacers and
at least two read pairs per edge. It generated 787 supplemental components
containing 4,040 copied rows (+119% array views, +56% rows) from 230 parents;
2,741 nodes remained isolated, and 246 degree>2 nodes appeared in 71 parents.
Under corrected strict scoring, CI changed from 15/9/14/20 to 20/4/14/20 and
CCF from 17/8/15/19 to 19/6/15/19; supported counts remained 957 and 1,047.
The two caller sets overlap, and these gains are grouping views of existing
spacers. Synthetic controls show a shared spacer can join two loci and a
single read spanning three spacers creates a nonadjacent co-occurrence edge.
Consequently these components are not safe to claim as independent arrays.
One focused merged family split into a 16-spacer component matching 16 of 17
reference spacers, but that in-sample example does not resolve the general
ambiguity. No component output was promoted.

Across the 661 primary parents, 3,965 co-occurrence edge records passed the
two-pair cutoff. Exact/RC matching against the high-confidence reference set
placed 981 at one annotated location with adjacent order, four at one
location with nonadjacent order, 91 at multiple possible locations or with
conflicting possible order, and two at no shared annotated array despite both
endpoints being annotated. Another 2,887 had at least one endpoint absent
from this reference subset and are unassessed. The two no-shared-annotation
edges are not established false positives; incomplete annotations and shared
spacers limit locus attribution. This analysis collapses caller overlaps by
source location and does not claim independent observations from both callers.

An independent synthetic string/graph control reproduced the same failure
modes: two disjoint three-node loci were recovered exactly, whereas a shared
spacer, a branch, or a repeated off-target co-occurrence merged distinct
loci; missing links fragmented a locus. A read containing three successive
spacers also yielded an edge between its first and last, which are not
adjacent. In 1,000 fixed-seed trials with two disjoint three-node loci and
four possible independent observations per junction, reducing the per-
observation junction probability from 1.0 to 0.5 reduced exact component
recall from 1.000 to 0.478. This is a conditional synthetic graph test, not
an estimate of biological detection performance.

The simulated read set lacks independent held-out genomes, strain mixtures,
varied abundance and true-negative calibration. Future trials need those,
full-span sequence accuracy, one-to-one locus assignment, correct order and
unsupported-join rates, separate unknown calls, and read-level provenance.

During this investigation the disposable scanner sources and outputs are in
`/tmp/my-mcaat-read-rescue.5ZhbQi/`, `/tmp/mcaat_all_links.cpp`,
`/tmp/my_mcaat_component_synthetic.py`, and `/tmp/mcaat-component-experiment/`.
These are scratch artifacts, not a
released or validated API. SHA-256 of the read-rescue source is
`914dc7a18ad2c5a83b8f7a470e92366e753c2eee1ebff1660b5bb837275d97e5`;
the allele-seeded scanner source is
`000479617b1b41a4fef06944f2cb881b591278f1d858322f33a3d8bc14f06773`;
the read-link scanner source is
`83a273794d532fb6465beef5c8cc972466e97d2068d2abd97461f05d31fc8c75`.
The synthetic-control script is
`8e5b69998357aff0fafa36ba526e964dd552ccc66d25bf7d8a2deafe17d65f31`.

### Graph cutoff and operational integrity

On the same retained graph, a four-thread default-cutoff control and a run
with the start-node multiplicity cutoff reduced from 20 to 10 generated 598
arrays / 6,711 rows and 609 arrays / 6,060 rows, respectively. Using the
high-confidence strict allele-aware audit, CI changed from 14/8/16/20 and
932 supported spacers to 14/9/15/20 and 928; CCF changed from 16/7/17/19
and 1,021 to 16/7/18/18 and 1,050. One CI partial became a miss and one CCF
split/merged call became partial. These are exploratory paired settings, not
causal effect estimates: dynamic start collection and per-thread visited
state also change enumerated cycles across runs. The new cutoff is not
promoted on this evidence.

An identical-settings replay of the default-cutoff control confirms that
this variability is material. With the same retained graph, image, four
threads and memory limit, the first and second runs produced 133,638 versus
136,314 cycles and 598 arrays / 6,711 rows versus 603 / 6,824. After
canonicalizing cycle rotation and reverse complement, they had 38,176 and
38,023 distinct cycles, of which 34,695 were shared (Jaccard 0.836). The
outputs had 5,640 and 5,750 distinct reverse-complement-canonical spacer
strings, of which 3,299 were shared (Jaccard 0.408). These are exact-string
comparisons, not false-positive or biological-loss rates: different phase or
boundaries can represent related evidence. The stronger change downstream is
consistent with sensitivity to input ordering, but does not establish its
mechanism. Comparisons to the singleton and cutoff ablations are consequently
single-trial observations, not reproducible effect estimates.
The singleton run had 91,072 distinct canonical cycles and shared 25,806 with
the first control (Jaccard 0.249), a larger observed change than between the
two identical-setting controls; this still cannot separate scheduling from
the parameter effect.

Strict allele-aware scoring also changed materially on the control replay:
CI moved from 14/8/16/20 to 17/9/16/16 (clean/split/partial/miss), with
932 versus 1,208 supported spacer occurrences; CCF moved from 16/7/17/19 to
19/8/18/14, with 1,021 versus 1,329. The two reference sets overlap, and
these are caller-based sequence matches rather than independent truth. Even
the aggregate performance is therefore unstable across these two nominally
identical runs, so small single-run differences should not be treated as a
validated improvement. Each caller had eight status improvements and one
regression. Most of the support increase comes from one region reported by
both callers: its CI count rose from 32 to 255 (223 of the 276 increase), and
its CCF count from 0 to 254 (254 of 308). This is one overlapping region, not
two independent improvements.

The audit also found that `low_abundance=true` does not prevent the initial
unconditional deletion of multiplicity-one edges. That interacts with tip
pruning, short repeats, branching thresholds, bounded enumeration, start
saturation, thread-dependent visited state, repeat grouping, spacer
clustering and boundary trimming. Cycle-only post-processing cannot recover
evidence lost earlier. A singleton-preserving graph run is an isolated
ablation, not evidence that preserving every singleton is safe: sequencing
errors can increase graph complexity and false candidates.

The isolated default-low-abundance ablation retained singleton edges before
tip pruning while leaving all other graph logic unchanged. On the same graph,
four threads and the same memory limit, its post-pruning graph retained
579,970,131 edges versus 297,839,652 in the control; start nodes increased
from 14,460 to 32,599 and retained cycles from 133,638 to 271,283. It emitted
1,376 arrays / 7,830 rows versus the control's 598 / 6,711, and took 991.4 s
versus 197.0 s end to end for graph replay. This exceeds the earlier 2x
runtime bound and more than doubles the number of candidate arrays. As with
the cutoff experiment, schedule-dependent enumeration is a confounder; these
single wall-clock runs also shared the machine with other analyses, so their
timings are not a controlled performance benchmark.

The ablation has 545 two-spacer and 278 three-spacer blocks, compared with
236 and 92 in the control. An exact/RC-normalized output comparison retained
396 of the control's 5,640 unique spacers, with 5,244 no longer exact and
6,655 additional candidate sequences. This measures output instability, not
necessarily biological loss, because boundary shifts can change strings.
For example, 2,642 of the 5,244 control-only exact strings have an exact
containment relation to an ablation output, and at least 1,602 have a seeded
ungapped match at 90% identity and 90% reciprocal span. These tests overlap,
are not one-to-one, and the seed can miss matches; they illustrate
representation churn without establishing correctness.
There were no malformed blocks, invalid bases, or sub-20-nt spacers, but the
ablation had two within-array RC-canonical duplicates versus one in the
control. The scratch build predates the later output-error fixes, which do not
change successful graph-stage behavior.

Under strict allele-aware scoring against the matched single-run control, CI
Bona-fide changed from 14/8/16/20 and 932 supported spacers to 1/3/5/49 and
33, while CCF Evidence-4 changed from 16/7/17/19 and 1,021 to 2/1/6/50 and
31. Both numbered ablation output files were included. For four representative
CCF clean-to-miss loci, an exhaustive strict scan of every ablation spacer
found zero matching spacers anywhere (0/37, 0/56, 0/180 and 0/24). Thus those
examples are not simply lost at the scorer's candidate-index or repeat gate;
altered boundaries or phase remain possible mechanisms. A permissive exact-
containment scan finds 32/37, 24/56, 99/180 and 24/24 of those sequences,
respectively, so the strict collapse must not be described as wholesale loss
of the underlying sequence. Nearly all of those containment hits are embedded
inside longer output spacers, supporting an overextension hypothesis without
proving locus identity. Only 0/32, 0/24, 1/99 and 1/24 of those containment
hits occurred in a repeat-compatible family, suggesting family assignment is
also involved. These are
single-trial, caller-based measurements, not general failure rates. A paired
status audit found 31 CI and 34 CCF regressions alongside two and one
improvements, respectively; overlapping callers are not independent loci. The
variant fails the runtime, output-growth and strict-score gates and was not
promoted.

The permissive allele-aware scorer also worsened: CI 22/16/14/6 and 1,439
supported spacers became 4/16/30/8 and 1,183; CCF 24/14/18/3 and 1,604 became
8/12/34/5 and 1,346. This score still permits exact containment at shorter
reciprocal span and should be interpreted alongside the strict result.

Separate from detection accuracy, a missing or unreadable cycles file and
failed array/cycle output writes could previously appear successful. The
postprocessor and cycle writer now propagate open, read, write and close
errors. A successful replay preserves its timestamp-normalized output, while
missing inputs, a directory used as input, an obstructed output path and a
full output device return errors. Experimental sidecar markers now carry
through file splits. The full image containing those code changes was
`sha256:3a2d428d3374e24df3a54f7b0a44d07b3fa43ec429dc44f337fd5fbac828a676`;
its frozen-cycle primary and sidecar normalized hashes remain the values
recorded above. CLI tests include a tiny real-graph cycle-file open failure.
The multi-file sidecar marker and a cycle-writer write/close failure were not
fault-injected. Remaining integrity risks include non-atomic replacement
of existing numbered output files and headerless cycle records that cannot
detect all DNA-only corruption or truncation; those require separate format
and publication designs, not silent compatibility-breaking changes.

## Implementation audit and new council — 2026-11-06

The new council split the work between graph traversal, post-processing, and
independent validation. Experiments are isolated under
`/tmp/my-mcaat-council5.hzz29U`; production changes preserve the prior dirty
worktree. The immediate plan is to establish a sorted-start control, compare
proved graph-search corrections on that control, and test circular
normalization together with both repeat boundaries rather than changing the
rotation alone. Array scores are accompanied by exact full-spacer matches
whose repeat is compatible in the same output row and orientation.

Fresh primary-source review found that the original bounded-cycle algorithm
has published incompleteness counterexamples:
[Bauernöppel and Sack, 2025](https://arxiv.org/abs/2512.08392) and their
[2026 bounded-scope algorithm](https://arxiv.org/abs/2607.14745).
The current MCAAT implementation differs from that pseudocode and passes the
five-node published counterexample, so the papers do not by themselves prove
this implementation loses those cycles. Independent exhaustive tests are
required before adopting a replacement.

The implementation audit found separate, reproducible mistakes: a node-only
visited set in the depth-limited prefilter suppresses a shorter return route
after a longer route; target closure at exactly the maximum length is checked
too late; relaxation temporarily narrows 64-bit node IDs to `int`; and a
cycle-closing neighbor remains on the enumeration stack and can be emitted
again after another branch. The last error can inflate support and consume
the 500-cycle cap. Search correctness and downstream detection effects are
tested separately below.

### Promoted experimental boundary mode

The paired representation correction is now available as
`--boundary-mode two-sided` or `boundary-mode=two-sided` in a settings file.
`legacy` remains the default. A mode marker is added only to two-sided output;
the existing supplemental repeat-view switch remains separate. The mode:

1. validates and removes the exact 22-base serialized closure;
2. counts circular-core 23-mers and resolves anchor ties lexically;
3. estimates the repeat on both sides of the anchor, using the existing 60%
   agreement threshold on distinct rotated cores and at least two distinct
   cores, with row-local repeat compatibility within two edits;
4. leaves at least 20 spacer bases and caps the combined repeat at 55 bases;
5. retains an exact shared-boundary fallback for rows that do not support the
   majority, instead of forcing the proposed repeat onto every row;
6. disables subsequent legacy shortening and tail correction, since the
   boundaries have already been inferred;
7. selects a row's repeat from observations of its chosen exact spacer, rather
   than borrowing a repeat supported by other members of a transitive spacer
   component.

The crucial error was coupled: the old rotation includes the redundant
closure, while its tail validator requires the repeat prefix before an anchor
to resemble the repeat suffix after the anchor. Those are generally different
parts of the repeat. Correcting only the rotation, as in the earlier rejected
experiment, leaves the incompatible boundary rule in place.

Frozen input: `ArraySearch/my_mcaat/cycles.txt`, SHA-256
`b1bcb5cf08fcdda3f41b511aeb169c68a8a18304049d6093c1ac1da897ce2a57`.
The 58 CI Bona-fide and 59 CCF Evidence-4 reference rows overlap; results are
caller-based and must not be added as independent biological discoveries.
Strict scoring uses the existing allele-aware comparator at 90% spacer
identity and reciprocal span, without its historical containment override.

| Same frozen input | CI clean/split/partial/miss | CI supported occurrences | CCF clean/split/partial/miss | CCF supported occurrences |
|---|---|---:|---|---:|
| Legacy | 15/9/14/20 | 957 | 17/8/15/19 | 1,047 |
| Final two-sided | 22/15/14/7 | 1,437 | 24/14/17/4 | 1,606 |

The final mode emits 680 arrays / 5,384 spacer rows versus the frozen legacy
661 / 7,212. More matches therefore do not come from simply appending many
alternative outputs. The independent exact-boundary audit requires an exact
full spacer and compatible repeat from the same output row, with consistent
orientation and an observed repeat allele in that reference. It deduplicates
canonical spacer sequences within each reference before summing:

| Exact association audit | Legacy | Final two-sided |
|---|---:|---:|
| CI full spacers with same-row repeat support | 655 / 1,519 | 1,119 / 1,519 |
| CCF full spacers with same-row repeat support | 762 / 1,703 | 1,568 / 1,703 |
| CI macro exact recall | 0.4114 | 0.6781 |
| CCF macro exact recall | 0.4729 | 0.7959 |

Of the final paired matches, all 1,119 CI and 1,567 CCF also match an exact
observed reference repeat; one CCF repeat is only approximately compatible.
This verifies family association, not spacer adjacency, order, or unique
locus identity. Every one of the 5,384 complete output row sequences
`repeat + spacer` equals an input circular core up to rotation/reverse
complement. Legacy has 2,532 of 7,212 rows without that complete-cycle match;
these are not automatically false biological calls, because its shortening
changes total sequence length. The invariant demonstrates sequence provenance
independently of caller labels.

Two independently generated control cycle inputs yielded exactly the same
strict status vectors as the final frozen result; supported occurrences were
1,437/1,436 CI and 1,605/1,605 CCF. Exact associations were 1,120/1,119 CI
and 1,569/1,567 CCF. The sorted-control, corrected-reachability, and corrected-
closure graph inputs also gave the same strict status vectors, with
1,436/1,437/1,437 CI and 1,604/1,606/1,606 CCF supported occurrences. Their
exact association counts were 1,120/1,119/1,119 CI and 1,568 throughout CCF.
The apparent clean-to-miss graph regression at `NC_005230_1` under legacy
boundaries disappears: all three corrected-boundary outputs recover its 45
distinct exact spacers together in a 45-row array.

The independent synthetic audit used four predetermined random seeds and
phase, strand, order and duplication perturbations. Each seed yielded
384/480 exact generating repeat/spacer pairs in the final approach, versus
268–276/480 for legacy. These repeated synthetic comparisons are not 480
independent biological loci. The remaining 96 failures involve divergent
repeat branches sharing an anchor. Random cycles, duplicate-only circles,
single-circle substitution errors and the tested low-complexity backgrounds
produced no candidates. Tandem-copy backgrounds and unrelated origins sharing
a 23-mer can still produce arrays; cycle sequence alone cannot resolve their
biological origin.

Promotion is deliberately opt-in. In the frozen paired audit, two overlapping
references in each caller change clean to split/merged even though their 19
and 11 spacers remain; two split calls become partial and some small partial
families become misses. Callers also disagree on a one-base boundary in one
large locus. Thus the large aggregate gains do not satisfy the earlier
absolute no-regression gate for a default replacement. Further work needs
independent biological truth, small-family vote calibration, divergent-repeat
branch separation, and explicit locus linkage. No threshold was tuned to a
particular reference. Full duplication/phase invariance and calibrated
specificity are not established.

Production legacy replay retains normalized SHA-256
`8c905ac91bb301c813003d2b3eb92959eb2d56fc27653560067e0b8a102be6c6`.
Production two-sided replay matches the tested scratch output exactly after
removing the timestamp and new mode marker, SHA-256
`169a40a4c8d56390371ce844338cc995df661633912119c02e21396004bf3bcd`.
One shared-host harness replay measured 7.80 s / 155,760 KiB maximum RSS for
two-sided versus 7.52 s / 180,912 KiB for legacy; these single runs do not
establish a speed improvement.

### Traversal fixes, rejected alternatives, and integrity

The bounded reachability prefilter now uses BFS with reusable thread-local
storage, which preserves shorter routes and accepts a closure at the maximum
length. Relaxation retains 64-bit node IDs. Enumeration consumes each closure
neighbor once and logs both cycle-cap and step-cap saturation. Start buckets
are sorted after parallel collection. Two four-thread sorted controls produced
identical cycle SHA-256
`5ef06c02af66329f44bbb1b37b9cb645b715bf91aa88151bb497c7e5884661db`.
This is reproducibility at fixed effective team size, not thread-count
independence.

On the retained 634,337,564-node graph, sorted legacy reachability admitted
14,460 starts; corrected reachability admitted 18,121. Before closure repair,
the sorted control contained 8,434 extra identical cycle records. With both
fixes, 108,926 records were all exact-unique; 61 starts reached the 500-cycle
cap and none reached the ten-million-step cap. Its cycle SHA-256 is
`904f31137041503c8d49eda14d16b850a1db1377b61f3ed9570729caab0e2901`.
Corrected graph output with two-sided boundaries is 664 arrays / 5,317 rows,
with the same strict status vectors as the frozen benchmark above.

The council also tested, and did not promote:

- Unanimous two-sided extension: some exact gains, but strict scores improved
  on one input and regressed on others. Existing repeat variants stop exact
  unanimous extension prematurely. Distinct-core majority votes addressed
  that mechanism without reference-specific thresholds.
- Global canonical-cycle deduplication: choosing one orientation split
  repeat-family evidence and reduced exact and strict recovery.
- Refreshing graph tips after singleton removal: removed all remaining tips
  and left final legacy arrays byte-identical, but changed 167 cycle records
  in each direction, all within two capped starts. There were 91 old-only and
  118 new-only canonical cycles. The single run took 224.6 s versus 189.3 s
  for the closure control on a shared host; no speed benefit was established.
- Fully independent starts, with cross-start visited suppression removed:
  completed in 540.1 s at a sampled 1.93 GiB, versus 189.3 s for the control.
  It produced 1,066,519 exact-unique records and 56,752 distinct circular/RC
  cycles, versus 108,926 and 36,903. It added 19,861 canonical cycles and lost
  12, all losses attributable to two starts reaching the 500-record cap.
  Corrected boundaries yielded 669 arrays / 5,483 rows. Exact reference
  associations changed CI 1,119 to 1,120 and CCF remained 1,568: one actual
  canonical spacer gained and one lost, with no new recovered reference or
  best-array gain. Strict status vectors were unchanged, with CI supported
  occurrences 1,437 to 1,438 and CCF unchanged at 1,606. That tradeoff does not
  justify the 2.85x replay cost.

An additional implementation audit removed a validation-error path that
offered recursive deletion of the output directory and read an uninitialized
answer on stdin EOF. Invalid settings now fail without deleting existing
output. This is separate from the remaining non-atomic array-file replacement.

Validation: 31,500 bounded-search oracle comparisons; 9,401 comparisons of
the production enumerator body against exhaustive simple-cycle search;
targeted duplicate/cap and high-node-ID counterexamples that fail on the old
logic; CLI exact-pair, phase, allele, input-format, negative, invalid-mode and
settings-override tests; and input/output fault and invalid-settings file
preservation checks. The built image is `my_mcaat:council5-final`; Python CLI
tests ran on the host because the builder has no Python interpreter. The
two C++ CTest targets passed in the build container. Existing compiler warnings
remain in vendored and unrelated code. Traversal allocations are not hard
bounded by `--ram`, and some graph counters still use 32-bit integers.

All experimental reports, scorecards and helper sources are consolidated under
the scratch root named above; only source, regression tests and this existing
log/README were added or updated in the repository. No new GitHub release was
created in this investigation.
