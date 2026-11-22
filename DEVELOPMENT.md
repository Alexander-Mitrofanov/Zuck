# Zuck development

## Supported-call output profiles — 14 November 2026

The native pipeline now exposes an additional selected-call layer while
preserving its inclusive primary rows. `call-profile=standard` selects
26–50 nt spacers with the existing observable exact-context support rule;
`context` omits the profile length bounds and `off` disables the new output.
Selection uses observed read evidence, never benchmark references or source
genomes. The ordinary confidence annotation and this selected view share
one raw-read evidence pass, including when the annotation table is disabled.

The four transaction-managed outputs are `call_decisions.tsv`,
`supported_entries.tsv`, `supported_arrays.tsv`, and
`supported_spacers.fasta`. Selection is per oriented row with original IDs.
The canonical spacer inventory accepts a sequence if any row passes, without
licensing other repeat pairings. Group membership retains singletons and
explicitly labels them. Missing read evidence remains unknown; unselected
rows and all original grouping sidecars remain available. Existing linked
or database subgroup memberships must be reconciled with selected original
row IDs before selected-array evaluation.

A fresh 35,132,388-read native run and a saved-cycle replay reproduce the
earlier diagnostic selection on this development community:

| Output | Canonical spacers | Exact union matches /2,053 | Exact union F1 | Strict clean CI /58; CCF /59 |
|---|---:|---:|---:|---:|
| Inclusive Zuck | 3,496 | 1,589 | 57.27% | 36; 36 |
| Supported standard | 1,757 | 1,586 | 83.25% | 31; 30 |
| Frozen CRASS comparison | 1,754 | 1,437 | 75.49% | 28; 28 |

The supported view selects 3,290 oriented rows in 133 original groups, with
54 explicit singleton groups. Evaluation credits only the 50 of 99 native
linked views whose entire original row membership survives; no partial path
is inferred. Strict union F1 is 85.46%, versus 84.79% for CRASS. The 1,739
demoted spacer types remain in the inclusive inventory: 1,304 are present in
source genomes (including three reference matches), and 435 lack an exact
source occurrence. Source presence does not establish a CRISPR role.

These are annotation-concordance results on one reused development community,
not held-out biological accuracy. The full run takes 601.3 seconds with
approximately 3.984 GiB largest-process peak RSS; one timing is descriptive,
not a speedup estimate. CRASS is not rerun in this experiment. The unchanged
original tables and native supported outputs agree between replay and full
inference; original numbered text differs only in its timestamp.

A separate deterministic raw-read stress fixture has 96 planted spacer types
at mixed lengths. Inclusive and `context` recover 80 types with 88 predictions
(F1 86.96%); `standard` recovers 56 with 64 predictions (70.00%), demoting 24
real 22- or 55-nt types. This exposes the length-profile tradeoff rather than
establishing generalization. The fixture contains random arrays and tandem
repeat controls, not representative biological communities.

All 34 CTest suites and four independent evaluator-integrity tests pass.
Coverage includes profile validation and precedence, row identity, length
endpoints, unknown evidence, duplicate/strand handling, singleton retention,
input protection, rollback, stale-file removal, empty output and preservation
of inclusive calls. Commands, frozen policy, source/binary identities, native
outputs, source-presence losses and scores are retained under
`ArraySearch/zuck_supported_calls_20261114/` in the enclosing workspace.

## Independent ±3 nt boundary candidates — 13 November 2026

The database-assisted pipeline and optional candidate diagnostics now expand
every exact pre-collapse oriented pair into all 49 combinations of independent
repeat-start/repeat-end shifts. A single boundary moves by at most three
bases; spacer length can change by six when both ends move. All bases remain
on the observed circular sequence. All combinations are exported when
candidate evidence is requested; discovery/output length constraints apply
to promotion. Duplicate seeds, orientations and source assemblies cannot
inflate corroboration.

The selected `catalogue` policy uses exact quality-qualified linear context
and known repeat/spacer membership in one eligible source array. It excludes
the focal spacer and circular core from its corroborating evidence and
requires a source-weighted score of at least 0.5. It takes the union of
unambiguous winners from nested ±1, ±2 and ±3 selection windows, then requires
at least two canonical spacer types per exact-repeat group. Existing calls
and linkage/database views remain intact. Alternative policies retain all
catalogue cuts, one full-window winner, or anchored novel spacers. Diagnostic
collection does not enable promotion with database read validation disabled.

The frozen, target-assembly-excluded development experiment evaluates
3,117,723 combinations from 63,627 seeds. There are 442,260 distinct full
context patterns and 134,070 eligible read-supported oriented pairs. The
comparison records 67 policy/window configurations, including controls.
The selected policy has the highest combined exact-reference F1 among those
tested; broader read-only and repeat-only rules inflate the candidate set.
The independent evaluation checks every combination against both exact
spacer panels after inference. Benchmark labels and source genomes are never
inputs to native selection.

| Metric | Previous database-assisted default | Selected ±3 policy |
|---|---:|---:|
| Exact CI spacers | 1,145 / 1,517 (75.478%) | 1,150 / 1,517 (75.808%) |
| Exact CCF spacers | 1,592 / 1,703 (93.482%) | 1,604 / 1,703 (94.187%) |
| Exact union spacers | 1,615 / 2,053 (78.665%) | 1,630 / 2,053 (79.396%) |
| Exact CI reference F1 | 44.388% | 44.410% |
| Exact CCF reference F1 | 59.570% | 59.795% |
| Exact union reference F1 | 56.716% | 57.043% |
| Strict CI / CCF matches | 1,442 / 1,614 | 1,443 / 1,615 |
| Strict CI / CCF reference F1 | 55.902% / 60.393% | 55.725% / 60.205% |
| Strict clean CI / CCF arrays | 37 / 37 | 37 / 37 |
| Primary groups / oriented rows | 681 / 6,605 | 687 / 6,710 |
| Distinct canonical spacers | 3,642 | 3,662 |

All 20 added spacer types and 16 added literal repeat+spacer types occur
exactly in the source genomes. No previous oriented row or source-supported
sequence is removed. Exact source presence does not validate a biological
boundary. Strict F1 declines slightly in CI/CCF despite increased strict
recall: the selected objective is exact boundary recovery with conservation
of existing calls, not improvement of every metric.

A fresh 35,132,388-read run reproduces the 163,622 original cycles and the
selected prototype's exact pair set. The final release is additionally
checked by complete-table replay, following a diagnostic-only guard fix.
Regression coverage is 32 passing CTest suites, including planted spacers
that need three-base extension or shortening at both ends, all 49 coordinates,
strand covariance, literal-context equivalence, self-corroboration rejection,
duplicate-source invariance, unrelated-source isolation and result integrity.

The full commands, immutable executables, source snapshots, per-variant
recall labels, policy tables, resource measurements and independent database
audit are preserved in `ArraySearch/zuck_boundary_20261113/`. See its
`REPORT.md` and `final_benchmark.json` for final measurements and limitations.
The initial screening uses raw spacer unions; the ranking/final prototypes
also enforce group-size constraints. Selection windows in that screen use
scores from the full ±3 universe; they are not independent native radius
runs. This is one repeatedly used development community, not held-out
validation, and reference F1 is concordance rather than biological precision.

## Code structure and review — 13 November 2026

The production entry point is `src/main.cpp`. It runs graph construction,
loading, cycle search, post-processing and cleanup through `core/pipeline.h`.
Debug builds use this same pipeline. The former genome-specific debug main
and its unused, incomplete filters/placeholders have been removed.

| Area | Responsibility |
|---|---|
| `core/cli.cpp`, `core/settings.h`, `core/option_values.h` | Resolve settings-file defaults and CLI overrides, validate complete values, preserve exact read paths, and protect inputs before writing. |
| `core/sdbg_build.cpp`, `core/graph_artifacts.h` | Adapt reads to MEGAHIT and define which generated files cleanup owns. Temporary aliases accommodate MEGAHIT's whitespace-delimited library format. |
| `crispr_array/cycle_finder.*`, `bounded_cycle_search.h`, `search_partitions.h` | Bounded graph traversal with deterministic logical partitions. |
| `crispr_array/post_processing.h`, `post_processing_helpers.h` | Named stages for cycle loading, boundary extraction, repeat clustering, spacer grouping, optional views and output. |
| `crispr_array/read_evidence.h` | Shared exact-pattern matching and lazy quality-span indexing for database, element and linkage evidence. Near-threshold spans use direct summation to preserve decisions. |
| `crispr_array/array_database*`, `database_groups.h`, `candidate_evidence.h`, `element_confidence.h`, `read_linked_groups.h` | Separate catalogue, candidate, element and adjacency evidence; scoring does not silently alter accepted calls. |
| `crispr_array/boundary_candidates.h` | Enumerate both junction shifts before collapse; share identical read contexts, exclude self-corroboration, select catalogue-supported alternatives, and write complete candidate provenance. |
| `io/result_bundle.h`, `io/run_manifest.h` | Stage result tables, publish the manifest last, and recover from rename errors. |
| `scripts/build_array_database.py` | Build verified, bounded catalogue projections and protect output/manifest paths. |

The remaining `spacer_ordering`, `core/main_run_and_debug`, older `io/reads`
and `protospacer_detection` helpers are experimental library code. They are
not invoked by the production CRISPR pipeline. Their edge-case regressions
check engineering behavior; they do not establish biological ordering or
protospacer validity. The obsolete `--benchmark` option is rejected; use the
external scoring workflow or `bench/run_benchmarks.sh` for graph timing.

`zuck_core` compiles production implementation files once for the executable
and C++ regression tests. Python tests cover the real CLI, graph construction,
replay, output preservation and shell runner. `ENABLE_DEBUG` enables optional
JSON diagnostic helpers while keeping the same entry point; it is distinct
from the standard `CMAKE_BUILD_TYPE=Debug` compiler configuration.

The council reviews first-party Zuck code and the two enclosing workspace
scripts (`tools/download_genomes.py`, `tools/run_crisprcasfinder.py`). Vendored
libraries, external caller copies, original reads and historical benchmark
artifacts are preserved. Pre-existing uncommitted edits are part of the
starting baseline, not changes attributed to this review.

The completed council validation passes 30 integrated CTest suites, six
offline workspace workflow tests, algorithm sanitizer/oracle checks and a
DEBUG JSON smoke test. Six full-corpus replays cover all boundary modes and
the optional view/filter/rescue settings; scientific tables and stage counts
match the saved baseline in every case. Single replay observations were
faster in all six cases, but shared-host timing is not a controlled speedup
estimate.

A fresh run over all 35,132,388 original reads completed in 676.289 seconds
with 4,177,280 KiB peak RSS. Its 163,622 cycles and all nine scientific TSVs
are byte-identical to the prior run; all stage counts match, retaining 681
groups and 6,605 rows. Numbered text output differs only in its timestamp.
This verifies compatibility on the existing development community, not new
biological validation. The full run froze the release build before late
DEBUG/experimental-helper cleanup; the final release rebuild passes the full
regression suite and six replay comparisons. A final indentation-only cleanup
has no semantic effect.

The separate patch, source snapshots, commands, checksums and council verdicts
are recorded in the enclosing workspace at
`ArraySearch/zuck_quality_council_20261113/`. The comparison policy for
human-readable output explicitly removes only its generation timestamp.

## Reproducibility and candidate evidence — 13 November 2026

The code/tool review identified three reproducible issues: search histories
changed with worker assignment, repeat/spacer dimension settings had no
effect, and singleton edges were removed despite `low-abundance=true`.
It also showed that transitive shared-23-mer collapse could discard one of
two planted internal-SNP spacer alleles, with representative choice changing
when only redundant traversal counts changed.

The first implementation package separates logical search partitions from
execution workers, enforces optional output constraints and fixes singleton
policy. Eight fixed partitions preserve the prior eight-worker development
search. Visited state resets at every logical partition, including when
OpenMP creates fewer workers than requested. The partition count remains an
inference parameter: this change fixes worker dependence, not exhaustive
cycle enumeration. Enumeration parallelism is limited to the partition count;
graph construction and start-node scanning retain their requested teams.
Singleton retention is opt-in; the new false default
matches the previous unconditional removal. Explicit old true settings now
retain singletons and must be changed to false for historical behavior.

Length bounds default to zero because the old nominal values were unused.
Explicit bounds filter both de novo and database rows before generating
row-indexed views. They constrain output, not candidate generation. This
preserves the established unbounded-by-settings default while making the
settings' actual effect inspectable in the manifest.

The optional candidate ledger scores every exact oriented de novo pair after
boundary extraction and before repeat-family/spacer collapse. It retains
traversal/core provenance separately from quality-qualified read evidence.
It does not use scores to choose representatives. This is the measurement
foundation for the review's proposed read-supported allele selection and
error-model experiments; multiplicity-dependent allele selection is still
present in the accepted-call pipeline. The ledger also cannot recover
candidates lost before boundary extraction or prove a biological boundary.

Regression coverage increases from 21 to 25 CTest suites. New checks include
288 production-search comparisons against serial partition histories under
normal, dynamic and thread-limited OpenMP teams; a real SDBG fixture across
1/2/4/8 workers; singleton-edge recovery; all four length limits; CLI/settings
precedence; database/confidence/linkage row alignment; and candidate bundle
rollback, missingness and stale-output removal. In the planted internal-SNP
fixture, all nine exact candidates receive read evidence while the eight
accepted calls remain unchanged. Doubling traversal observations leaves
read scores unchanged. The real graph fixture yields 32 cycle records with
singleton removal and 36 with retention, including four records containing
the planted rare spacer only in the retention arm.

A fresh full-read run with candidate evidence enabled processed 35,132,388
reads in 716.264 seconds (11.94 minutes), with 3.984 GiB peak RSS. Its 163,622
cycle records and all eight existing TSV tables are byte-identical to the
archived confidence-enabled database run. Every prior stage count matches;
the output retains 681 primary groups and 6,605 rows. Consequently the
previous reference-recovery measurements are preserved; this package does
not establish a biological accuracy improvement. These are single-run
resource measurements, not an isolated estimate of candidate-scoring cost.

The ledger contains 63,627 exact candidates. Of these, 4,273 have at least
two quality-qualified sequence units and two context start offsets, including
530 exact oriented pairs absent from the primary output. These may represent
alternate cuts, rejected groups or discarded alleles. Read context alone
does not justify promoting them. A future selection experiment must compare
such candidates using explicit error/boundary models and independent recovery
and conservation gates rather than treating all 530 as new discoveries.

The review probes and fresh implementation experiments are preserved in the
enclosing workspace under `ArraySearch/zuck_code_review_20261113/` and
`ArraySearch/zuck_foundation_20261113/`. The latter includes the pre-change
source snapshot and executable, regression logs, full-run command and
machine-readable compatibility checks. Existing benchmark outputs are
preserved.

## Element and boundary confidence — 13 November 2026

The selected implementation adds `element_confidence.tsv` by default and
preserves the database release's complete calls. `--confidence-scores false`
disables annotation. Scores distinguish exact read presence of repeats and
spacers, literal repeat–spacer–next-repeat context, repeat-family support,
catalogue evidence, and ambiguity among emitted boundary decompositions.
They are bounded evidence indices, not calibrated correctness probabilities.
Unknown evidence is explicit. Complete read sequence units and start offsets
are capped with saturation flags; they are not proven independent molecules.

A new council considered quality/error models, duplicate dependence,
catalogue confidence/frequency, leave-one-out corroboration, competing cuts,
graph branches/coverage, paired/long reads, resampling, endpoint models, and
learned calibration. Five selection proposals were tested sequentially on
the frozen target-assembly-excluded database benchmark. The predeclared
default gate required preservation of reference recovery, clean arrays, and
source-supported spacer and literal-pair evidence.

| Policy | Exact CI / CCF | Strict clean CI / CCF | Lost source-present spacer / pair types |
|---|---:|---:|---:|
| Complete calls + confidence, selected | 1145 / 1592 | 37 / 37 | 0 / 0 |
| Unique catalogue boundary | 1090 / 1524 | 24 / 27 | 123 / 97 |
| Full read-context filter | 1144 / 1590 | 36 / 34 | 895 / 226 |
| Catalogue corroboration excluding focal spacer/core | 1112 / 1544 | 29 / 31 | 95 / 72 |
| Other-core repeat-family preference | 1145 / 1590 | 36 / 36 | 111 / 88 |
| Read + catalogue + family selection | 1089 / 1521 | 24 / 26 | 1047 / 341 |

The read-context view improves exact reference F1 substantially,
0.443884/0.595697 to 0.596921/0.791242, while losing supported sequences and
some reference matches and clean arrays. This is an explicit reference-concordance/coverage
tradeoff, not an across-the-board improvement. The complete calls remain the
default; the tested filtered view and its decision ledger are retained in the
experiment artifacts. Reference-panel F1 is not biological precision.

Read-presence indices outperform the tested catalogue/family alternatives for
their corresponding sequence-presence endpoints: AUROC 0.986486 for repeats,
0.986207 for spacers, and 0.976789 for literal repeat+spacer pairs. Evaluation
deduplicates canonical sequence types. Source presence does not establish a
CRISPR role or boundary. Catalogue-derived boundary ranking remains a useful
diagnostic, but its unique preference disagreed with exact same-array
reference membership in 33/56 eligible CI and 48/69 CCF contested classes.
Both panels are incomplete computational annotations; these are disagreement
counts, not independently established error rates.

All 6,605 rows receive an annotation, including explicit unresolved statuses.
The seven previous primary/evidence/grouping tables remain byte-identical.
The measured confidence-off replay took 85.086 seconds; annotation took
146.902 seconds, a 61.816-second increase on this one run. Both preserve
3,642 spacer types, exact recovery 1,145/1,592, strict recovery 1,442/1,614,
and strict-clean recovery 37/37. Tests include 21 CTest suites, 21 independent
native synthetic controls, seven rank-metric controls and 14 selection controls.
The fresh full-read run completed in 631.677 seconds (10.53 minutes), using
3.984 GiB peak RSS and 3,591.67 CPU seconds. All 130 final audit checks passed;
full-run cycles, confidence table, preserved result tables and benchmark
panels match replay. The existing prepared catalogue was reused. These are
single-run measurements, not a repeated performance estimate.
The full benchmark, audits, source identities, council notes and
all rejected proposals are recorded in the enclosing workspace's
`ArraySearch/zuck_confidence_20261113/REPORT.md`.

## Database guidance — 13 November 2026

The opt-in `--array-database` workflow now uses CRISPRidentify array confidence,
distinct same-array repeat/spacer matches, and exact quality-checked linear
read context. The default database policy uses Bona-fide sources, preserves
all de novo evidence, and emits maximal catalogue-coherent content views.
Possible-array guidance and novel-spacer recovery are explicit optional arms.
See `README.md` for preparation and inference commands.

Two councils designed and reviewed the implementation. The first integration
improved recovery but exposed circular-versus-linear adjacency errors. The
second implementation validates the full repeat, spacer, and 16 bases of the
next repeat before counting source-array matches. It also separates database
confidence tiers and adds source-coherent grouping without adding sequences
solely because they exist in the catalogue.

The primary development benchmark excludes all 45 target assemblies and their
GenBank counterparts from the database. Exact spacer matches improve from
1,121 to 1,145 of 1,517 CI types and from 1,570 to 1,592 of 1,703 CCF types.
Strict matching improves from 1,439/1,611 to 1,442/1,614; strict clean arrays
improve from 36/36 to 37/37. All 5,430 original oriented rows and the existing
read-linkage tables are retained. All 146 added spacer types and 150 added
literal repeat–spacer types have exact source-genome support. The 29 database
content views reference existing rows and do not establish genomic loci.

Exact F1 decreases from 0.447237 to 0.443884 in CI and from 0.603962 to
0.595697 in CCF; strict F1 also decreases. These are supported recovery and
grouping gains with a reference-concordance tradeoff. Source presence and
read adjacency do not establish CRISPR identity, and the reference panels
are incomplete. Related genomes remain after accession exclusion.

The fresh full run completed in 572.16 seconds (9.54 minutes), using eight
threads and 4 GiB requested graph RAM, at 3.98 GiB peak RSS. Its cycles,
primary calls, database evidence, grouping tables, and scores reproduce the
selected replay exactly. Database projection cost 164.80 seconds separately;
it is sample-dependent and requires cycles from a first pass. The native
full-run timing assumes the projection is already prepared. Single-run
timings do not support a robust speed comparison.

Nineteen CTest suites, ten database-view adapter tests, thirteen independent
read-context checks, and both full-catalogue maximal-view reconstruction audits
pass. The final full-run audit passes all 76 checks. The complete report and
all measured alternatives are retained in the workspace at
`ArraySearch/zuck_db_20261113/REPORT.md`.

## Validated de novo baseline — 12 November 2026

The follow-up council selected two-sided inference with additive terminal-variant
rescue and read-linked content views enabled, strand-equivalent content groups
coalesced, and short-spacer filtering disabled. Native replay on the same cycles
and 17,566,194 paired reads completes in 42.28 seconds, including the read scan.
The fresh full-read run completed in 540.55 seconds (9.01 minutes), at 3.98 GiB
peak RSS, with eight threads and 4 GiB requested graph RAM. It regenerated the
frozen cycles byte-for-byte and exactly reproduced replay primary/grouping
tables and scores. Fourteen CTest suites, ten adapter tests and all 121 final
audit checks pass. The active build matches the tested binary.

Compared with the prior validated default, exact recovery rises from 1,119 to
1,121 of 1,517 CI spacers and from 1,568 to 1,570 of 1,703 CCF spacers. Strict
recovery rises to 1,439/1,517 and 1,611/1,703. Strict clean-array recovery rises
from 25 to 36 of 58 CI arrays and from 27 to 36 of 59 CCF arrays. The 99 supported
views supplement all 610 original parent groups. Every prior oriented row is
retained, with nine added rows and six additional canonical spacer types.
Source-supported spacer types increase from 3,055 to 3,061 and source-supported
literal repeat-plus-spacer types from 3,950 to 3,957, with zero previous losses.
Two new literal pairs lack source support; source presence alone does not prove
CRISPR identity or array ownership.

Grouping alone reaches 37/37 strict clean arrays. Terminal rescue introduces
known competing endpoint alleles; the integrated default conservatively retains
that ambiguity and reaches 36/36 while improving spacer recovery. The grouping
maxima and recovery maxima must not be combined into an unmeasured claim. A
full-context uniqueness refinement improved one secondary CCF best-group match
but no primary recovery or clean-array endpoint, so the selected implementation
was retained. Independent community validation remains necessary.

Detailed results, improvement lists and rejected experiments are in the workspace
report `ArraySearch/zuck_recovery_round2_20261112/REPORT.md`. The prior full default
completed in 497.75 seconds at 3.98 GiB peak RSS; its unchanged artifacts remain in
`ArraySearch/zuck_improvement_20261112/`.

## Starting evidence: 12 November 2026

The latest completed benchmark used 17,566,194 paired reads from one simulated
community containing 45 assemblies and 99 replicons. It compared full raw-read
workflows and shared-assembly callers, with explicit replay timings. Frozen
CRISPRidentify Bona-fide (CI) and CRISPRCasFinder Evidence-4 (CCF) annotations
provide separate, incomplete computational reference panels.

| Boundary mode | CI exact spacer recall | CCF exact spacer recall | CI clean groups | CCF clean groups | Distinct spacer types |
|---|---:|---:|---:|---:|---:|
| legacy | 666/1,517 | 779/1,703 | 17/58 | 18/59 | 5,836 |
| two-sided | 1,119/1,517 | 1,568/1,703 | 25/58 | 27/59 | 3,490 |
| branch-aware | 1,119/1,517 | 1,568/1,703 | 25/58 | 27/59 | 3,545 |

The table describes the earlier boundary-mode comparison. Two-sided mode was
the strongest completed configuration under those measured
recovery and concordance endpoints. Strict spacer recall is 1,437/1,517 CI and
1,609/1,703 CCF for both two-sided and branch-aware. A clean group requires
at least 90% distinct reference-spacer recall and 80% output-group purity in
a repeat-compatible group. These groups are not validated genomic loci.

The legacy full run took 491.93 seconds; full run plus two-sided replay took
approximately 8.25 minutes. A direct two-sided full run was not part of that
measurement. Timing includes one fresh replicate on the development host.

## Current improvement priorities

1. Validate recovery and read-linkage rules on independent communities and
   sequencing-quality profiles before generalizing the development gains.
2. Recover candidates absent from retained cycles. Only five missed exact
   spacers in each panel existed among the pre-collapse candidate alleles;
   changing collapse alone cannot recover absent candidates.
3. Resolve terminal-allele and shared-spacer ambiguity using additional observed
   context. Do not hide known alternatives or enumerate arbitrary favorable
   subpaths from branched components to improve grouping scores.
4. Investigate unannotated read-supported repeat-spacer-repeat content. Broad
   rescue added 369 source-supported types but reduced reference F1 and is not
   enabled. Source presence alone does not prove CRISPR identity.
5. Review boundary discrepancies: exact-to-strict gaps remain 318 CI and 41 CCF
   spacer types. Separate reference disagreements from tool errors.
6. Inspect the 435 spacer types without exact source support and the two added
   literal pairs lacking source support.
7. Preserve recovery and same-row sequence evidence when changing grouping,
   duplicate handling, or low-support filtering.

Candidate changes must use the same frozen inputs, reference panels, and
scoring definitions. Record failures and rejected variants. A default should
be selected on explicit measured endpoints rather than a configuration name.
Changes selected on this community still need independent validation.

The default uses two-sided reconstruction, terminal rescue and read-linked
grouping enabled, and the short-spacer filter off.
The optional filter applies only to two-sided/branch-aware modes. Source
support revealed a real evidence loss despite unchanged high-confidence
reference recall, so incomplete-panel concordance alone did not determine
selection. Final coalescing validation preserved every oriented row and
the same recovery, concordance, same-row support, and clean-group endpoints.
There is no separate performance-mode flag: omitting `--boundary-mode`
selects the chosen default directly.

## Existing engineering constraints

Cycle replay is the fast path for controlled post-processing comparisons.
Use a fresh full-read run to validate end-to-end behavior after selecting a
candidate. Preserve the original benchmark inputs and output directories;
new experiments require new output paths and source fingerprints.

The branch-aware mode improves several predetermined synthetic mixtures but
has no additional exact reference recovery on this community. It can change
other supported output sequences, so equal panel recall does not imply equal
biological behavior. Alternative repeat views expose already accepted evidence
and must not be scored as new spacer discovery.

Logical partition count and visited history can change graph-search output;
execution worker count and effective OpenMP team size do not change those histories.
The graph RAM budget is not a hard traversal-memory limit. Stage observations
can contain duplicate traversals and must not be interpreted as read depth.
Groups do not establish unique loci, terminal repeats, or biological order.

The project retains source authorship and dependency notices. Detailed
pre-rename development records and rejected historical prototypes are kept
unaltered in the workspace provenance archive. The active documentation
describes current behavior and the present evaluation.
