# Zuck 7.0 spacer benchmark

Native adaptive Zuck had higher exact spacer F1 than both fixed Zuck 6.1 and CRASS 1.0.1 in all four low-depth comparisons. Tolerant F1 depended on the cohort and depth. CRASS led native Zuck at 3× in the new communities and at 1× in the earlier communities.

| Cohort | Coverage | Fixed 6.1 exact / tolerant | Native 7.0 exact / tolerant | CRASS exact / tolerant | 7.0 + CRASS exact / tolerant |
| --- | ---: | ---: | ---: | ---: | ---: |
| New communities | 1× | 52.19 / 54.28 | 54.87 / 58.45 | 47.90 / 53.36 | 55.64 / 59.14 |
| New communities | 3× | 75.38 / 78.44 | 76.33 / 80.19 | 72.36 / 83.75 | 76.54 / 83.33 |
| Earlier communities | 1× | 39.17 / 39.97 | 43.21 / 45.26 | 39.40 / 46.81 | 47.92 / 51.96 |
| Earlier communities | 3× | 68.63 / 71.84 | 71.06 / 76.11 | 64.95 / 75.99 | 72.48 / 79.90 |

Values are pooled F1 percentages. Each cell lists exact F1 first. The assisted column uses an external CRASS installation. [Counts, precision, recall and F1 by annotation panel](benchmark-7.0.tsv)

## Libraries and reference spacers

Three new five-genome communities were selected with no genome overlap with the 75 assemblies used in the earlier studies. For each community, a simulated 3× HiSeq paired-end library with 126-base reads supplied a nested 1× library. Each new-cohort row pools three libraries. The earlier cohort used three other five-genome communities and two draws per depth, giving six libraries per row. Each caller received the same reads. The 7.0 launcher inferred its route from FASTQ. The nominal coverage labels were used for evaluation.

The reference was the union of CRISPRidentify Bona-fide and CRISPRCasFinder Evidence 4 spacer sequences from the source assemblies. When the annotators disagreed, both distinct sequences entered the union. Identical sequences present in both were counted once per library. The linked table also reports scores against each panel separately. These reference annotations can miss real spacers and disagree on their endpoints.

## Scoring

Within each library, spacer sequences were deduplicated across strands. An **exact** match required identity over the entire sequence, allowing reverse complementation. A **tolerant** match required at least 90% ungapped identity and an overlap covering at least 90% of each spacer. Matching was one-to-one, so two predictions could earn at most one match to the same reference spacer. Counts were pooled across libraries before calculating F1 as twice the matched count divided by the sum of predicted and reference counts.

On the new 3× libraries, native Zuck made 281 predictions with 249 tolerant matches. CRASS made 257 predictions with 250 matches. The number of additional unmatched calls explains most of the tolerant F1 gap. Native Zuck had 237 exact matches and CRASS had 216.

One 8×, one 20× and one mixed-coverage control selected the full 6.1 route and reproduced its spacer inventory. These controls check fallback behavior. The comparison with CRASS at low depth applies to the simulated communities and sequencing model studied here.
