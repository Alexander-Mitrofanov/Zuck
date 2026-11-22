# Reading Zuck results

Each numbered text file describes a sequence-content group with a repeat
consensus and rows of repeat/spacer pairs. The matching row in `arrays.tsv`
summarizes the group. `entries.tsv` retains every reported pair exactly.
`entry_index` is presentation order, not inferred biological spacer order.

The optional [reference annotation stage](database-and-protospacer-matches.md)
adds exact repeat/spacer catalogue matches and candidate viral protospacers.
The web view highlights every matched sequence and links each source record;
additional matches remain available through pagination and the complete export.

Two-sided and branch-aware modes consolidate strand-equivalent group boundaries
without removing or reorienting any observed row. A group can therefore contain
both row orientations and multiple rows for one canonical spacer. These are
preserved observations, not new spacer types or extra independent support.

Repeat variants are reconstructed from graph cycles. Neither the consensus
nor a row proves that the sequence came from a particular genomic locus or
repeat copy. Repeated rows and traversal observations are not read-depth
measurements. Supplemental alternative views carry an explicit kind and
parent identifier and reuse existing sequence evidence.

Read-linked groups are additional views in `read_linked_groups.tsv`, not
replacements for the primary groups. Each membership references an exact
one-based parent `entry_index`, including all original row orientations for
each selected spacer. `read_linked_edges.tsv` records undirected spacer
connections supported by exact repeat/flank junctions. A full repeat and at
least 16 bases from each adjacent spacer must match; all available flank
bases up to the full spacer length are checked. A flank must identify one
canonical spacer within its parent group.

Each matched span requires Phred quality at least 10 at every base and an
expected-error sum of at most one (`sum(10^(-Q/10))`). This applies across the
full repeat and all verified flank bases. The rule admits otherwise exact
junctions with moderate individual qualities while bounding the aggregate
error estimate; it does not provide a calibrated probability for a group.

Support is counted across distinct full read sequences (one input file), or
distinct ordered pairs of full read sequences (two files). Exact duplicate
sequences count once, and each edge also needs two junction-start offsets.
These units are not independent molecules. Complete unbranched components
with at least three spacer types become views; branched and circular
components are not cut into arbitrary subsets. Parents remain available,
including evidence outside a supported view. Neither a path nor a view
establishes unique locus ownership or biological spacer order.

Full read runs use their input reads for this step by default. Cycles/graph
replays require `--grouping-reads`; without reads the manifest records that
linkage was not applied. FASTA records are counted as lacking quality and do
not support edges. FASTQ input errors abort replacement of the result bundle.

Use `run_manifest.json` to identify the compiled source, effective settings,
cycle input, and result files. A replay describes only the supplied cycles;
its current graph-search defaults do not establish how those cycles were
generated. Preserve the external benchmark version and reference strata
alongside any downstream scores.

Compare modes on the same cycle input for a controlled post-processing
experiment, then measure a fresh read-to-result run for end-to-end cost.
Incomplete computational annotations measure concordance and recovery,
not biological precision or specificity.

## Graph-seeded read recruitment

`--read-recruitment true` enables an additional read pass after graph-family
inference. It requires two-sided or branch-aware boundaries and FASTQ input
(`--input-files`, or `--grouping-reads` during cycle replay). It is disabled
by default. It searches reads using graph-derived repeat seeds, verifies
repeat copies with at most one substitution, and retains the exact observed
spacer interior. This bounded implementation handles repeat lengths 20–55
and spacer lengths 26–50 bases, subject to any narrower explicit output
length limits. It can append a missing repeat/spacer pair to one unambiguous
existing family; it does not create seedless families or establish a locus.

Three evidence kinds are recorded: `full_pair` observes two complete repeats;
`full_partial` observes a complete upstream repeat and at least 16 downstream
repeat bases at a physical read end; `graph_junction` validates an existing
graph proposal with at least 16 repeat bases on either side of a complete
spacer. All available bases of each verified partial repeat are checked.
None of these modes fills in an unobserved spacer. Each verified span needs
Phred quality at least 10 and expected errors at most one.

Evidence is pooled by the exact oriented repeat/spacer pair using set union,
with at least two distinct full read or paired-read sequence units and two
spacer-start offsets (measured 16 bases upstream of the spacer). Reverse
complements and mate swaps do not create new units. Counts saturate at 32;
unit identities use two independently seeded 64-bit hashes. Every possible
owning family must have at least two independently supported spacer types.
These units are not independent molecules or calibrated probabilities.

`read_recruitment_candidates.tsv` retains raw repeat/spacer/repeat observations
and directed pair aggregates. Raw triplets are provenance only; aggregates
list their source candidate IDs and observed downstream repeat alleles.
An aggregate never licenses its reverse pair implicitly. Boundary conflicts
are inherited from every source. Alternative cuts of an existing circular
unit remain evidence only, so this pass cannot repair an already incorrect
graph boundary. Ambiguous family ownership blocks new rows. Resource limits
discard affected recruitment evidence rather than selecting the first rows
encountered; original graph output remains available.

`read_recruitment.tsv` binds eligible aggregates to exact array and entry IDs.
Original IDs are preserved and new rows append deterministically. Under the
standard or context profile, original rows may be selected by literal context
or by recruitment evidence that passes the row-promotion policy below. Newly
observed pairs remain inclusive and cannot enter selected output. `call_decisions.tsv` identifies the channel and
candidate ID explicitly. Existing context measurements retain their original
meaning. The standard 26–50-base spacer rule still applies. The manifest
records the contract, counts, settings and source fingerprint. Disabling
recruitment restores the original selection rule and removes stale sidecars.

The tandem screen introduced in schema 3 separates observation from automatic
selection and adds a competing
tandem-repeat explanation. It tests each physical raw triplet using a motif
learned from a fully observed repeat, trying the left repeat and then the right
repeat in the reversed frame. A graph junction with neither repeat fully seen
can instead use its explicitly labeled graph repeat proposal. Every observed
flank base and the whole spacer must fit the same continuing phase with
similarity strictly greater than 0.82. The spacer allows edit distance; the
observed flanks use Hamming distance. The two models are a short periodic motif
and a complete repeat copy. Repeat complexity is diagnostic, not a blanket veto:
a repetitive repeat with diverse spacers can still supply eligible evidence.

An aggregate needs at least two sequence units and two oriented start offsets
from raw sources that this tandem model does not explain. Both counts are set
unions, not subtraction from total support. A physical tandem explanation is
shared by both directed pair emissions. The ledger records actual verified
flank lengths, the chosen training side and orientation, motif, fit scores,
and separate non-tandem support. This is a conditional sequence explanation,
not a probability that an array is biological or nonbiological. Repeat/flank
indels and tandem patterns outside the bounded models can remain undetected.

Schema 4 adds a conservative policy for promoting recruitment evidence. The
candidate discovery, tandem screen, and error calculation above are unchanged.
Promotion requires an exact oriented repeat/spacer row already present in the
same graph family before recruitment. For every declared owner, another
candidate must have eligible read evidence and an actual proof binding to an
original row in that owner. Its repeat must have equal length and differ by at
most one substitution in either strand. Its spacer must differ from the focal
spacer by at least 18% of the longer length, using the smaller Levenshtein
distance to either strand. Context candidates are fixed before promotion;
they do not have to pass this new policy themselves. The first qualifying
candidate ID, then original row index, is recorded deterministically.

`read_recruitment_original_rows.tsv` snapshots every original row, including
duplicates. `read_recruitment_context.tsv` records a witness or an explicit
missing witness for every owner of every bound candidate. The proof table
separates candidate evidence eligibility from final row eligibility. New pairs
have status `novel_boundary_unresolved`; original pairs with otherwise eligible
evidence but missing diverse context have `insufficient_spacer_diversity`.
Original literal-context calls remain available. Novel observations cannot
supply family-context repeat alleles, alter prior read-linked groups, or become
selected through a supplemental view.

This policy can miss real novel spacers and rare alleles absent from the
original graph output. These remain visible in inclusive output; observation
does not resolve their biological boundary. The policy is experimental, and
retrospective improvements do not establish superiority on new communities.

The error screen introduced in schema 2 remains. An observed exact
allele with an equally long, one-substitution spacer neighbor under the same
oriented repeat is screened for error ambiguity. A second read pass counts
every quality-qualified window consisting of the repeat's last 16 bases,
a complete spacer-length window, and 16 arbitrary following bases, in either
orientation. This deliberately broad opportunity count includes other spacer
interiors, duplicate records and possible alternative errors. The lowest
spacer-base Phred score across those windows supplies a conservative error
rate under a calibrated-quality assumption.

The screen compares `min(1, choose(N,K) * p^K)` with `0.05/M`, using the
opportunity count `N`, capped distinct **non-tandem** support `K`, worst quality-derived
error rate `p`, and `M = 3 * sum(spacer lengths)` over all admissible directed
pairs before either ambiguity screen. Tandem-ambiguous pairs retain their
exposure metadata but do not undergo the error calculation. It computes in log
space and considers every eligible
one-substitution neighbor without iterative pruning. The ledger records these
inputs, the chosen neighbor ID, status and selection eligibility. Missing or
resource-limited opportunity evidence blocks recruitment promotion.

This is a conditional screen, not a calibrated biological probability:
sequence units need not be independent molecules, sequencing errors may be
correlated, and broad opportunity windows can make the test conservative.
It can leave a true rare allele unresolved. Ambiguous alleles remain in the
inclusive output and evidence ledger. Original rows keep their original
literal-context selection; newly appended ambiguous rows cannot bypass the
screen through the later literal-context pass. No allele is merged into its
more common neighbor.
