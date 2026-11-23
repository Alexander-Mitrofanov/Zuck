# Adaptive calling in Zuck 7.0

Run the adaptive launcher from the extracted release directory. For a source build, use the second command from the repository root after building with CMake.

```sh
python3 zuck_depth.py adaptive --input-files R1.fastq.gz R2.fastq.gz --output-folder result
```

```sh
python3 standalone/zuck_depth.py adaptive --zuck build/Zuck --scanner build/zuck-kmer-spectrum \
  --input-files R1.fastq.gz R2.fastq.gz --output-folder result
```

The launcher estimates depth from the reads. A uniform estimate below 5× enables `--read-first-sparse-junctions true`. Intermediate, high, mixed and uncertain predictions use the full 6.1 logic. Both routes enable clipped-repeat rescue and family admission. The direct `Zuck` executable remains available for fixed-parameter runs. [Coverage model and legacy routing](COVERAGE.md)

At low depth a complete 26–50 nt spacer on a single read can qualify if oriented repeat segments flank both boundaries and other fragments contain two distinct spacer contexts with compatible geometry. A short repeat segment at a physical read end requires an exact hook and a complete exact opposite repeat. Longer segments can tolerate a substitution away from the junction. The caller checks read qualities and withholds proposals with competing boundaries or unsupported near-neighbour variants.

| Output | Interpretation |
| --- | --- |
| `read_first_spacers.fasta` | Native spacer inventory for this library |
| `read_first_entries.tsv` | Call origins, including `sparse_junction` |
| `read_first_sparse.tsv` | Per-proposal decisions on the low-depth route |
| `coverage_routing.json` | Estimated regime, depth, selected route, effective switches and exit status |
| `run_manifest.json` | Native version, settings and source fingerprint |

## Optional CRASS-assisted inventory

With CRASS 1.0.1 installed separately, run the following from the extracted bundle after a successful native adaptive run. Keep the original input reads available.

```sh
python3 zuck_sparse_assisted.py --native result --crass /path/to/crass --output-folder assisted
```

The merger invokes CRASS only for the low-depth route. It retains native spacer types and adds CRASS spacers with distinct interiors. It treats exact duplicates and short terminal shifts as already represented. The combined inventory is `assisted/assisted_spacers.fasta`, and `crass_decisions.tsv` records the treatment of each CRASS sequence. In a source checkout the script lives in `standalone/`. [Native and assisted benchmark results](https://github.com/Alexander-Mitrofanov/Zuck/blob/v7.0/docs/benchmark-7.0.md)
