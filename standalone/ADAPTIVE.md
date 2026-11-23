# Adaptive Zuck research caller

Run the native adaptive caller with the same estimated coverage model and input files used for a regular preflight.

```sh
python3 standalone/zuck_depth.py adaptive --zuck build/Zuck \
  --settings standalone/v6.1.conf --scanner standalone/zuck-kmer-spectrum \
  --input-files R1.fastq.gz R2.fastq.gz --output-folder result
```

The low-depth route (uniform estimate below 5×) enables `--read-first-sparse-junctions true` and retains the complete 6.1 family admission and shallow rescue settings. Intermediate, high, mixed and uncertain predictions keep sparse recovery off. `estimate` and the original `run` mode remain available. The new option requires a native binary built from this research revision. The main spacer FASTA is `read_first_spacers.fasta`. The native entry table labels additions `sparse_junction`. The low route writes a separate `read_first_sparse.tsv` decision ledger. The wrapper records predicted depth, route, native command and return code in `coverage_routing.json`. [Depth estimator and scanner build](COVERAGE.md)

Sparse recovery scans both strands for two compatible repeat segments flanking a complete 26–50 base read subsequence. Other fragments must show two dissimilar spacer contexts of compatible lengths. A short repeat hook at a physical read end requires exact sequence and stronger length evidence. Longer segments allow a mismatch away from the junction and use FASTQ qualities to check spacer and boundary reliability. The pass appends qualifying new sequence types after standard selection. Already-reported spacers, competing boundaries, ambiguous families, tandem explanations and unsupported singleton variants are withheld. The count is a distinct observed spacer type per library. All sequences in the FASTA have a native source independent of CRASS.
