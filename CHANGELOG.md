# MY_MCAAT changelog

## v1.1.0 — 2026-11-06

- Add experimental, opt-in `--boundary-mode two-sided` inference over circular
  cycle cores. Preserve observed repeat/spacer pairing, infer both repeat
  boundaries, and avoid applying legacy shortening a second time. Legacy
  boundary inference remains the default.
- On the same frozen cycle input with strict allele-aware scoring, clean
  reference calls improve from 15 to 22 (CI) and 17 to 24 (CCF); misses fall
  from 20 to 7 and 19 to 4. Exact full-spacer matches with same-row repeat
  support improve from 655/1,519 to 1,119/1,519 and from 762/1,703 to
  1,568/1,703. These overlapping caller references are not independent truth;
  small-family and split-call regressions remain. See `DEVELOPMENT.md` for
  scoring definitions, controls, synthetic tests, and limitations.
- Fix bounded graph reachability, inclusive maximum-length closure, 64-bit
  node handling, and repeated cycle-closure emission. Report cycle/step cap
  saturation and sort collected starts for fixed-team reproducibility.
- Add lossless compact cycle serialization, reducing one frozen cycle file
  from approximately 422 MB to 23 MB while preserving legacy replay output.
- Correct MEGAHIT RAM-unit conversion and automatic memory sizing. Traversal
  memory is still not hard-capped by `--ram`.
- Fail explicitly on unreadable cycle input and output open/write/close
  failures. Invalid settings preserve existing output instead of prompting
  for directory deletion. Array-file replacement is still non-atomic.
- Retain optional experimental alternative repeat views, which expose
  already-accepted allele evidence rather than recovering new spacers.
- Add five regression-test suites, including 31,500 bounded-search oracle
  comparisons, 9,401 production-enumerator comparisons, boundary/provenance
  tests, and CLI fault/preservation tests. Consolidate council decisions,
  literature review, rejected experiments, and benchmark caveats in the
  existing development log.

This is a source release. Benchmark inputs and scratch experiment/scoring
helpers are not bundled. No general runtime speedup or calibrated biological
specificity is claimed.

## Upstream changelog snapshot

This file records pre-fork MCAAT development notes and is retained only for
provenance. The current `my_mcaat` council iterations and benchmark results are
recorded in `DEVELOPMENT.md`.

## [Unreleased] - 0.5.1
- Fixed crash in spacer ordering when no reads are found (guard added).
- Improved parallel scaling:
  - Replaced frequent `#pragma omp critical` usage with per-thread buffers and serial merges.
  - Introduced a lock-free visited bitmap (atomic 64-bit words) to remove synchronization hot-spots.
- Reduced allocator contention and reused per-thread containers to lower memory churn under heavy parallelism.
- Build fixes: added missing includes and small portability fixes so CMake build succeeds on target platforms.
- Build: CMake configure and full build completed successfully (targets `mcaat` and `runTests` built).


*Notes:* these changes focus on improving scalability and preventing serialization on large, memory-bound graphs.
