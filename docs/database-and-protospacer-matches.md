# Database matches and candidate protospacers

`scripts/annotate_results.py` adds source-linked evidence to a completed Zuck
result. It reads **every repeat consensus, repeat variant and spacer** in
`arrays.tsv` and `entries.tsv`, deduplicates sequence searches, and retains the
binding back to every occurrence through its sequence identifier. It does not
change array reconstruction, boundaries, confidence scores or grouping.

## Methods reviewed

* [CRISPRloci (2021)](https://pmc.ncbi.nlm.nih.gov/articles/PMC8265192/) includes
  viral input searched against a spacer collection. Its
  [CRISPRVirusSpacers implementation](https://github.com/BackofenLab/CRISPRloci/blob/main/CRISPRVirusSpacers/CRISPRVirusSpacers.py)
  uses `blastn-short`, word size 7 and an E-value cutoff of `1e-7`. That search
  direction differs from finding viral targets for Zuck's reconstructed spacers.
* [CRISPRTarget](https://pmc.ncbi.nlm.nih.gov/articles/3737339/) searches supplied
  spacers against nucleotide databases and supports inspection of target flanks
  and PAMs. It is a relevant alternative; PAM scoring is outside this integration.
* CRISPRloci v2 already uses
  [CRISPRspacer](https://github.com/Alexander-Mitrofanov/CRISPRspacer) `viral-search`
  for the required spacer-to-reference direction. Zuck reuses this standalone
  engine and its native evidence replay command. The tested version is **0.3.1**,
  commit **ddea3f5da9d5ff6a33802046bdb0e0a762029d23** (MIT).

The adopted search uses local `blastn-short`, word size 7, both strands, no dust
or soft masking, ungapped alignments, reward 1, penalty -3 and retrieval E-value
1000. The acceptance filter requires **100% query coverage and at most two
substitutions** (configurable to zero or one in the standalone annotator).
These are the existing CRISPRloci v2 settings, not the older CRISPRloci settings.
The reference record count is passed as the BLAST target limit. The engine's
limits fail explicitly instead of publishing truncated evidence.

BLAST retrieval is heuristic, even at zero substitutions: terminal mismatches
can be clipped and subsequently rejected. No accepted hit does not establish
absence. A candidate match does not establish infection, host range or active
targeting. PAM compatibility and overlap with CRISPR arrays are not evaluated.
Spacers outside 18–80 nt or containing ambiguous DNA are explicitly ineligible.

## Catalogue matching

Use the immutable CRISPRloci snapshot, including `catalog.sqlite3`, its
`exact-search.sqlite3`, and `shards/`. The registry digest and shard inventory
must match the completed search index. The index nominates candidate shards
for full-length sequences on both strands; all nominated shards are verified
by checksum and decompressed with a size bound. Actual sequence equality is
required after nomination. A hash match alone never establishes evidence.

Repeats are checked against both array consensuses and individual repeats.
Spacers are checked against observed spacers. All matching source records are
retained with category, certainty score, source tool, array, accession and
coordinates when supplied. These are sequence annotations: the inference
thresholds used by `--array-database` do not suppress them. Duplicate source
occurrences are not independent confidence votes. Ambiguous query sequences
are marked ineligible for exact catalogue matching.

Alternatively, `--array-database arrays.tsv` searches a supplied native TSV
projection. Its restricted scope is recorded explicitly; the sample-conditioned
projection is not a substitute for a full catalogue. Source links point to
NCBI nucleotide intervals, or assembly accessions when that is the only
source identity supplied by the projection.

## Run

Install the pinned CRISPRspacer release in a dedicated environment and install
local NCBI BLAST+. Provision a validated immutable viral reference using
CRISPRspacer's `viral-index` / `validate-viral-reference` operator workflow.
Ordinary Zuck jobs perform no downloads or reference updates.

```bash
python3 scripts/annotate_results.py --results result \
  --catalogue /references/crisprloci-snapshot \
  --viral-reference /references/viral-refseq \
  --viral-sha256 <sha256-of-reference.json> \
  --crisprspacer /opt/crisprspacer/bin/crisprspacer \
  --blastn /usr/bin/blastn --threads 4 --timeout 1800
```

Either search can run independently. An omitted reference produces
`not_configured`, never a negative search result. A supplied but invalid
reference fails the annotation. Existing output is never overwritten; use a
fresh `--output` directory to run an alternative annotation.

The atomic `annotations/` output contains:

* `queries.tsv`: sequence IDs, kind, search states, complete match counts and
  ineligibility reasons.
* `matches.jsonl`: every source match and candidate protospacer, with links and
  alignment evidence. Join `query_id` to `queries.tsv`, then kind/sequence to
  each native row or consensus.
* `manifest.json`: native input hashes, annotation hashes, reference identity,
  methods, tool version, BLAST identity/settings and limitations.
* `spacers.fasta` and `native-protospacers/`: submitted unique spacer sequences
  and the engine's original replayable evidence, when viral search runs.

The web download includes the first four flat files (FASTA when present).
Its JSONL contains all accepted evidence, while private native BLAST scratch
files are removed with the job's working directory. The web index verifies
annotation hashes, native input identity, sequence coverage and match counts
before publishing any result. Older results remain readable and explicitly
show that reference searches were not run.

Resource limits are 100,000 distinct annotation sequences and 1,000,000 total
matches; the reused CRISPRspacer engine additionally applies its own query,
raw-HSP and accepted-HSP limits. Exceeding a limit fails without completed
annotations. The worker also retains cancellation, disk, memory and time limits.

## Validation

```bash
python3 tests/test_annotations.py
ZUCK_CRISPRSPACER_TEST=/opt/crisprspacer/bin/crisprspacer \
  python3 tests/test_annotations.py
```

Controls cover independent repeat/spacer matches, reverse complements,
consensus-only matches, candidate hash false positives, corrupted shards,
complete source retention, native forward/reverse BLAST matches, substitutions
and rejected partial alignments. Web tests additionally cover row/match
pagination, authenticated source access, downloads, missing references and a
real queued native/BLAST job.
