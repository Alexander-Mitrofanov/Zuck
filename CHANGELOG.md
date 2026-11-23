# Zuck changelog

## Zuck 6.1 research candidate — 23 November 2026

- Added optional recruitment of complete spacers where a read ends inside a
  repeat. Independent spacer contexts establish the repeat family before
  these observations can add a call.
- Added family-specific screening of single-fragment read-only spacer
  variants when the graph and independent fragments already support a
  nearby sequence. The evidence ledger records the selection decision.
- Added `standalone/v6.1.conf` to run both changes with the 6.0 benchmark
  settings. Use the two switches independently for ablation studies.

## Zuck 6.0 — 23 November 2026

- Added `--read-first true` to discover repeat families in FASTQ reads,
  verify spacer junctions and combine read-supported spacers with the graph
  inventory in `read_first_spacers.fasta`. Companion TSVs record the source
  of each call and template support for each family.
- Added `--adaptive-coverage true` to measure repeated graph-edge
  multiplicities and adjust cycle discovery in sparsely supported libraries.
  The run manifest records the decision and the effective search settings.
- Added `--read-first-family-recut true` to correct a short graph repeat
  boundary when independent fragments support distinct, longer repeat
  families. The per-occurrence evidence table records the decision and
  direct counts of Q20 whole repeat–spacer–repeat fragments. All three new
  switches default to off in the native executable.
- Benchmarked the frozen release candidate against CRASS 1.0.1 on 30 held-out
  libraries. Zuck's exact spacer F1 was higher by 3.30–12.97 percentage
  points across the five coverage profiles. CRASS's tolerant F1 was higher
  by 1.26–2.90 points. The family correction changed zero inventories in
  those libraries and improved exact recovery from 46/223 to 199/223 in a
  separate targeted 80× related-repeat mixture. [Command, metrics and results](README.md#read-first-spacer-discovery)

## Standalone container artifact — 16 November 2026

- Build and test the Docker image as a portable, fully static executable in a
  `scratch` runtime, and add a BuildKit `standalone` export containing the exact
  tested binary, checksums, licensing, and attribution.
- Add an ELF-level standalone test that rejects dynamic loaders and shared
  library sections, then exercises the CLI with an empty executable search path.

## Experimental family repeat context — 16 November 2026

- Add a separate opt-in `--repeat-mixture-refinement true`, requiring
  `branch-aware` boundaries, for bounded nested family partitioning and
  five-base allele signatures. Defaults remain unchanged. This post hoc model
  repairs tested repeat mixtures but fails correlated spacer-terminal
  controls (40/40 to 0/40 exact calls); it is not a general accuracy upgrade.

- Add opt-in `--family-repeat-context true` to validate an exact upstream
  repeat–spacer against a compatible observed next-repeat allele in the same
  output family. Keep the default exact-same-allele rule unchanged.
- Preserve sequence-unit, offset, quality and missingness checks; record the
  selected context policy in the manifest. Discovery and inclusive rows are
  unchanged.
- Add controlled mixed-repeat, strand, wrong-family, quality, duplicate and CLI
  regression cases. These establish local behavior, not full-cohort superiority.

## Unreleased — 14 November 2026

- Add default `--call-profile standard` supported-call outputs using existing
  quality-qualified full-context evidence and an explicit 26–50 nt profile.
  `context` keeps any candidate length; `off` disables the selected view.
- Preserve all inclusive primary rows and original IDs. Export per-row
  decisions, selected oriented entries, family memberships with labelled
  singletons, and a canonical spacer FASTA. Keep missing evidence unknown.
- Share the confidence read scan with selection. `--confidence-scores false`
  disables its annotation table; use `--call-profile off` as well to disable
  this evidence pass when no other enabled feature needs it.
- Manage the four new files in the result transaction, including input-alias
  protection, rollback and removal of stale outputs when profiles are disabled.
- Pass 34 CTest suites and reproduce the supported profile in a fresh full-read
  run: development-set exact union F1 83.25%, with three fewer exact matches
  than the inclusive output. Preserve all original evidence. Mixed-length
  synthetic tests expose the standard profile's sensitivity cost; `context`
  remains available without the 26–50 nt bound.

## Unreleased — 13 November 2026

- Expand every pre-collapse repeat/spacer candidate into all 49 independent
  ±3 nt junction combinations for database-assisted inference and optional
  candidate diagnostics. Score literal read contexts and exact catalogue
  membership; retain source-array corroboration with the focal spacer/core
  excluded. Add nested-window catalogue selection, explicit experimental
  policies, coordinate/evidence sidecars and transactional cleanup.
- Cover both shorter and longer planted spacers, including three bases at
  each end, strand covariance, duplicate evidence, unrelated-source isolation,
  and diagnostic-mode isolation. All 32 CTest suites pass. See
  `DEVELOPMENT.md` and the workspace boundary experiment for measured results.

- Complete an iterative first-party quality council review. Split CLI handling
  and post-processing into named components; share exact matching and quality
  checks, use bounded edit distance and avoid unused optional-view work.
- Reject malformed CLI/settings values consistently and preserve paths with
  spaces. Protect input aliases, stage cycle output, and restrict cleanup to
  generated graph artifacts; preserve supplied graphs and unrelated files.
- Repair evaluation, ordering and experimental path edge cases. Use the same
  production entry point for debug/release; remove obsolete incomplete code
  and reject the previously ineffective production `--benchmark` option.
- Expand integrated coverage to 30 suites, including strict settings, input
  preservation, algorithm oracles, experimental paths and the benchmark runner.
  Add an optional DEBUG JSON suite and six offline workspace workflow tests.

- Separate logical cycle-search partitions from OpenMP workers. Keep eight
  fixed partitions by default to preserve the prior eight-worker search;
  expose `--traversal-partitions` as an inference setting.
- Make `low-abundance=true` retain singleton graph edges. Default to false
  to preserve the previous effective behavior; earlier true settings removed
  singletons unconditionally and now intentionally behave differently.
- Enforce optional repeat/spacer output length bounds from CLI and settings
  files before generating row-indexed sidecars. Zero disables a bound.
- Add opt-in pre-collapse candidate allele evidence with independent raw-read
  scoring, explicit primary membership and traversal/core provenance. Keep
  allele selection and accepted calls unchanged.
- Pass 25 CTest suites and a fresh 35.1-million-read run. Preserve the archived
  eight-worker baseline's 163,622 cycles, all eight existing TSV tables, 681
  groups and 6,605 rows exactly. Record 63,627 diagnostic candidates, including
  530 read-supported pairs absent from the primary output, without promoting
  them to calls.
- Add complete per-entry repeat, spacer, linear-context and boundary evidence
  annotations, including missingness, bounded counts, source provenance and
  explicit competing-cut ambiguity. Enable by default; expose
  `--confidence-scores false` for the original inference cost.
- Preserve all current calls after five confidence-guided filtering proposals
  failed the declared recovery/source-conservation gates. Retain the strongest
  F1 tradeoff as a separately documented experimental view.
- Extend native regression coverage to 21 CTest suites. Add independent
  planted-boundary, identical-context ambiguity and score-ranking controls.
- Add opt-in CRISPRidentify array-database guidance with bounded, confidence-
  weighted distinct repeat/spacer evidence and same-array co-occurrence.
- Prepare sample-specific catalogue projections by streaming and checksum-
  verifying source shards. Preserve confidence and provenance; deduplicate
  array signatures and support counterpart-accession exclusions.
- Require complete observed cycle sequences and quality-checked linear
  repeat–spacer–next-repeat context before database-supported discovery.
  Keep Bona-fide guidance as the default database policy; expose Possible
  evidence and novel-spacer recovery as explicit experimental options.
- Add maximal source-array content views with validated row membership,
  subset suppression, and bounded ambiguity. Preserve all original de novo
  rows and read-link views.
- Record database identity, effective policy, source signatures, distinct
  match counts, read support, and content-view provenance in the result bundle.
- The database integration added regression coverage to 19 CTest suites, including catalogue
  preparation, read-context corroboration, and database grouping.

## Unreleased — 12 November 2026

- Adopt the Zuck name for the executable, source identifiers, build/install
  targets, manifests, result headers, and active documentation.
- Preserve original authorship, license terms, and dependency notices.
- Complete a council-led evaluation against the frozen 12 November benchmark;
  follow-up recovery/grouping work passes fourteen CTest suites and ten
  grouping-adapter tests. The fresh full run matches replay exactly and passes
  all 121 final audit checks. See `DEVELOPMENT.md` for measured results and
  rejected alternatives.
- Make two-sided boundary reconstruction the default and retain explicit
  legacy/branch-aware overrides in the CLI and settings files.
- Add an opt-in short-spacer filter for modern boundary modes. Keep it off by
  default: its development ablation removes eight source-supported spacers,
  so higher incomplete-panel F1 does not establish a better default.
- Coalesce redundant strand-equivalent group boundaries while retaining every
  observed oriented repeat/spacer row. Groups remain unordered sequence
  content and do not establish genomic loci.
- Enable bounded terminal-variant recovery, adding observed equally supported
  one-base terminal alleles without replacing accepted rows.
- Enable read-linked content views using exact junctions, bounded expected
  errors, distinct sequence units and start diversity. Preserve parents and
  expose exact row membership and support in transactional TSV sidecars.
- Reject truncated gzip inputs and input/output aliases before replacing
  results; remove stale read-link sidecars when grouping is disabled.

## Existing capabilities

- Replay retained compact cycle evidence without rebuilding the graph.
- Infer both repeat boundaries with `two-sided`, with optional experimental
  `branch-aware` separation of divergent families sharing an anchor.
- Preserve observed repeat/spacer pairing and bounded repeat-family
  aggregation; optionally emit labelled alternative-repeat views.
- Apply bounded reachability checks, inclusive maximum cycle length,
  64-bit graph node handling, and duplicate cycle-closure suppression.
- Write deterministic array/entry tables and a versioned manifest containing
  source identity, settings, stage/filter counts, and checksums.
- Stage result bundles, roll back rename failures, reject unsafe output
  destinations, and preserve earlier results when parsing or staging fails.
- Validate graph search, boundary reconstruction, CLI behavior, and result
  integrity through the CTest suites.

Historical development notes and pre-rename documentation remain in the
workspace provenance archive rather than being presented as new releases.
