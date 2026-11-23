# Zuck 7.0 for Linux

Extract the Linux x86-64 release archive, enter the extracted directory and run the commands below. The caller and scanner are statically linked. Coverage estimation and adaptive calling use Python 3 and its standard library.

```sh
sha256sum -c SHA256SUMS
./Zuck --version
python3 zuck_depth.py adaptive --input-files reads_R1.fastq.gz reads_R2.fastq.gz --output-folder result
```

Use one FASTQ for single-end reads or two mate files for paired-end reads. Both plain and gzip-compressed FASTQ are supported. The benchmark preset uses four threads and a 4G native RAM setting. Choose a new or empty output directory.

The native spacer inventory is `result/read_first_spacers.fasta`. Check `result/coverage_routing.json` to see which calling mode was selected. [Adaptive calling and optional CRASS support](ADAPTIVE.md) · [Coverage estimation](COVERAGE.md) · [Benchmark results](https://github.com/Alexander-Mitrofanov/Zuck/blob/v7.0/docs/benchmark-7.0.md)
