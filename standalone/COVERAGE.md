# Coverage estimation in Zuck 7.0

The coverage preflight takes one FASTQ or a pair of Phred+33 FASTQ files, plain or gzip. It estimates base sequencing depth for a roughly uniform shotgun sample and flags strong evidence of mixed depth. It uses the input reads alone. An uncertain or mixed sample has a null `depth_x`. Mixed-model component weights refer to repeated 31-mer types. They cannot identify taxa, the fraction of a genome covered, or the depth at a particular CRISPR locus.

## Run the preflight

The Linux bundle includes the scanner, launcher and presets. Run the following from its extracted directory to print the estimated regime and depth.

```sh
python3 zuck_depth.py estimate --input-files reads_R1.fastq.gz reads_R2.fastq.gz
```

From a source build, use `python3 standalone/zuck_depth.py estimate --scanner build/zuck-kmer-spectrum --input-files reads_R1.fastq.gz reads_R2.fastq.gz`.

Use the `adaptive` command for Zuck 7.0 calling. The `run` command reproduces the earlier 6.1 shallow-only policy. The following table describes that legacy policy. [Zuck 7.0 calling modes](ADAPTIVE.md)

The launcher prints JSON. A native run also writes `result/coverage_routing.json` with the prediction, underlying spectrum, chosen arguments, explicit overrides and exit status. Use `--scanner` to select a separately built scanner.

| Prediction | Sequencing depth model | Zuck settings |
| --- | --- | --- |
| low | Uniform estimate below 5× | Enable clipped-repeat rescue and disable dense family admission |
| intermediate | Uniform estimate from 5× to below 12× | Full 6.1 preset |
| high | Uniform estimate at least 12× | Full 6.1 preset |
| mixed | Two separated depth components improve fit | Full 6.1 preset |
| uncertain | Insufficient repeated sequences, low quality, out-of-range depth or poor fit | Full 6.1 preset |

The routing thresholds describe this research configuration. To override either switch explicitly, supply `--read-first-shallow-rescue true|false` or `--read-first-family-admission true|false` to `run`. Use the exact v6.1 preset for comparison with the study.

## Interpretation

The scanner builds a bounded hash sample of distinct strand-canonical 31-mers and counts all observations of retained words. The fitted Poisson rate is corrected for eligible windows, base quality and read length. Very small inputs can yield an uncertain result. Genomes too rare to produce enough repeated words may be invisible. Copy-number structure, PCR duplicates, overlapping mates, sequencing bias and quality-score calibration can change the spectrum. Validation used simulated paired-end HiSeq short reads from three five-genome communities.
