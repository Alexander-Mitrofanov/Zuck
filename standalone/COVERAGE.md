# Coverage estimation and optional Zuck 6.1 routing

The coverage preflight takes one FASTQ or a pair of Phred+33 FASTQ files, plain or gzip. It estimates base sequencing depth for a roughly uniform shotgun sample and flags strong evidence of mixed depth. It uses the input reads alone. An uncertain or mixed sample has a null `depth_x`. Mixed-model component weights refer to repeated 31-mer types. They cannot identify taxa, the fraction of a genome covered, or the depth at a particular CRISPR locus.

## Run the preflight

In the private Linux x86-64 research bundle, `zuck-kmer-spectrum`, `Zuck`, `v6.1.conf` and the two Python modules are in the same directory. The launcher uses Python 3 with its standard library. Run these commands from the extracted bundle directory.

```sh
sha256sum -c SHA256SUMS
python3 zuck_depth.py estimate --input-files reads_R1.fastq.gz reads_R2.fastq.gz
python3 zuck_depth.py run --input-files reads_R1.fastq.gz reads_R2.fastq.gz --output-folder result
```

From a source checkout, build the companion scanner with a C++17 compiler and zlib development files. Run this from the Zuck source root.

```sh
c++ -O3 -std=c++17 standalone/zuck_kmer_spectrum.cpp -lz -o standalone/zuck-kmer-spectrum
python3 standalone/zuck_depth.py estimate --input-files reads_R1.fastq.gz reads_R2.fastq.gz
```

To estimate depth and run the existing Zuck 6.1 executable with the selected options, use an empty or new result directory.

```sh
python3 standalone/zuck_depth.py run --zuck build/Zuck --settings standalone/v6.1.conf \
  --input-files reads_R1.fastq.gz reads_R2.fastq.gz --output-folder result
```

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
