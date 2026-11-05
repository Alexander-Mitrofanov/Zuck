# my_mcaat development log

This is the single living record for council decisions, implementation
iterations, and benchmark outcomes. Disposable builds and runs stay outside
the repository under `/tmp`.

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
