# Benchmark

## Test tools

All tools support for random access were used to perform benchmark. Python 3.6.6 was used to test tools that was developed as Python package.

| Tools    | Version  | Language  | URL                                       |
|----------|----------|-----------|-------------------------------------------|
| samtools | v1.9     | C         | http://www.htslib.org/                    |
| pysam    | v0.15.3  | C, Python | https://github.com/pysam-developers/pysam |
| seqkit   | v0.11.0  | Go        | https://bioinf.shenwei.me/seqkit/         |
| pyfasta  | v0.5.2   | Python    | https://github.com/brentp/pyfasta         |
| pyfaidx  | v0.5.5.2 | Python    | https://github.com/mdshw5/pyfaidx         |

## Test data

We downloaded 15 genome FASTA files from NCBI FTP server with different size and sequence counts. These genomes were listed in table below. 

| Species name             | Common name      | Assembly Accession | Genome size (bp) | Sequence counts |
|--------------------------|------------------|--------------------|------------------|-----------------|
| Ambystoma mexicanum      | Axolotl          | GCA_002915635.1    | 32,393,605,577   | 125,724         |
| Pinus lambertiana        | Sugar pine       | GCA_001447015.2    | 27,602,651,501   | 4,253,096       |
| Picea glauca             | White spruce     | GCA_000411955.5    | 24,633,074,982   | 3,033,321       |
| Triticum aestivum        | Bread wheat      | GCA_900519105.1    | 14,547,261,565   | 22              |
| Triticum dicoccoides     | Wild wheat       | GCA_900184675.1    | 10,494,955,245   | 149,145         |
| Rana catesbeiana         | Bullfrog         | GCA_002284835.2    | 6,250,335,504    | 1,544,634       |
| Locusta migratoria       | Migratory locust | GCA_000516895.1    | 5,759,798,599    | 1,397,492       |
| Pisum sativum            | Pea              | GCA_003013575.1    | 4,275,931,177    | 5,449,423       |
| Homo sapiens             | Human            | GCF_000001405.38   | 3,257,302,968    | 593             |
| Zea mays                 | Maize            | GCF_000005005.2    | 2,134,513,431    | 266             |
| Gallus gallus            | Chicken          | GCF_000002315.5    | 1,065,348,650    | 463             |
| Sorghum bicolor          | Sorghum          | GCF_000003195.3    | 708,876,072      | 868             |
| Cicer arietinum          | Chickpea         | GCF_000331145.1    | 530,893,862      | 7,127           |
| Arabidopsis thaliana     | Thale cress      | GCF_000001735.4    | 119,300,826      | 6               |
| Saccharomyces cerevisiae | Yeast            | GCF_000146045.2    | 12,071,326       | 16              |

Download these FASTA files, and place in a directory (e.g. data), uncompress the file and retain the gzip compressed file.

## Test method

We used linux command ``/usr/bin/time -f "%e %M"`` to estimate the running time and peak memory of building index, random access and reverse complement for each tool. Each test was performed in triplicate. The mean was calculated as final results.

## Perform test

### Building index

```sh
./benchmark_build_index.sh 3 data/*.fa > benchmark_result_build_index.tsv
python make_ggplot2_matrix_for_index_access.py benchmark_result_build_index.tsv > build_index_matrix.tsv
```

### Random access to longest sequence

```sh
./benchmark_random_access.sh 3 data/*.fa > benchmark_result_random_access.tsv
python make_ggplot2_matrix_for_index_access.py benchmark_result_random_access.tsv > random_access_matrix.tsv
```

### Reverse complement longest sequence

```sh
./benchmark_reverse_complement.sh 3 data/*.fa > benchmark_result_reverse_complement.tsv
python make_ggplot2_matrix_for_reverse_complement.py benchmark_result_reverse_complement.tsv > reverse_complement_matrix.tsv
```

## Plot

Recommended you use RStudio to plot, replace the input file in Line 4 of Plots.R file with matrix file generated above, and then copy code to RStudio to plot.
