# Benchmark

## Test tools

All tools support for random access were used to perform benchmark. Python 3.6.6 was used to test tools that was developed as Python package.

| Tools    | Version  | Language  | URL                                         |
|----------|----------|-----------|---------------------------------------------|
| pyfastx  | v0.7.0   | C, Python | <http://github.com/lmdu/pyfastx>            |
| samtools | v1.9     | C         | <http://www.htslib.org/>                    |
| pysam    | v0.15.3  | C, Python | <https://github.com/pysam-developers/pysam> |
| seqkit   | v0.13.2  | Go        | <https://bioinf.shenwei.me/seqkit>          |
| bioperl  | v1.7.7   | Perl      | <https://bioperl.org>                       |
| pyfasta  | v0.5.2   | Python    | <https://github.com/brentp/pyfasta>         |
| biopython| v1.74    | Python    | <https://biopython.org>                     |
| pyfaidx  | v0.5.8   | Python    | <https://github.com/mdshw5/pyfaidx>         |

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

Download these FASTA files, and place in a directory (e.g. data/fastas), uncompress the file and retain the gzip compressed file.

We downloaded 5 FASTQ files from [NGDC database](https://bigd.big.ac.cn/). These FASTQ files were listed in table below.

| Accession | Total bases \(bp\) | Read counts  | Plain file size \(B\) | Plain file size \(GB\) | Gzip file size \(B\) | Gzip file size \(GB\) |
|-----------|--------------------|--------------|-----------------------|------------------------|----------------------|-----------------------|
| CRD000339 | 39,221,236,800     | 392,212,368  | 102,558,579,168       | 95\.52                 | 30,316,624,332       | 28\.23                |
| CRD000389 | 20,655,260,429     | 204,507,529  | 49,907,114,035        | 46\.48                 | 17,716,170,588       | 16\.50                |
| CRD000371 | 13,930,424,000     | 139,304,240  | 36,426,113,582        | 33\.92                 | 10,972,470,560       | 10\.22                |
| CRD000411 | 5,515,159,742      | 54,605,542   | 13,325,681,150        | 12\.41                 | 4,204,994,577        | 3\.92                 |
| CRD000365 | 1,597,691,427      | 15,818,727   | 3,860,387,755         | 3\.60                  | 1,211,374,918        | 1\.13                 |

Download these FASTQ files, and place in a directory (e.g. data/fastqs), uncompress the file and retain the gzip compressed file.

## Test method

We used linux command ``/usr/bin/time -f "%e %M"`` to estimate the running time and peak memory of building index, random access and sequence iteration for each tool. Each test was performed 3 times. The average value was calculated as final results.

## Perform test

### Building index

building index for fastas

```sh
./benchmark_fasta_build_index.sh 3 data/fastas/*.fa > benchmark_result_fasta_build_index.tsv
python3 make_fasta_ggplot2_matrix.py benchmark_result_fasta_build_index.tsv > matrix_fasta_build_index.tsv
```
building index for fastqs

```sh
./benchmark_fastq_build_index.sh 3 data/fastqs/*.fq > benchmark_result_fastq_build_index.tsv
python3 make_fastq_ggplot2_matrix.py benchmark_result_fastq_build_index.tsv > matrix_fastq_build_index.tsv
```

### Random access

randomly get 30% of sequence names from fastas with random seed of 123

```sh
python3 get_random_sequence_names_from_fasta.py 123 0.3 data/fastas/*.fa
```
random access to sequences from fastas

```sh
./benchmark_fasta_random_access.sh 3 data/fastas/*.fa > benchmark_result_fasta_random_access.tsv
python3 make_fasta_ggplot2_matrix.py benchmark_result_fasta_random_access.tsv > matrix_fasta_random_access.tsv
```

randomly generate one thousand subsequence region with length of 1 Kb and random seed 123

```sh
python3 generate_random_interval.py 123 data/fastas/*.fa
```

random access to subsequences from fastas

```sh
./benchmark_fasta_extract_subsequences.sh 3 data/fastas/*.fa > benchmark_result_fasta_extract_subsequences.tsv
python3 make_fasta_ggplot2_matrix.py benchmark_result_fasta_extract_subsequences.tsv > matrix_fasta_extract_subsequences.tsv
```

randomly get ten thousand names of reads from fastqs with random seed of 123

```sh
python3 get_random_read_names_from_fastq.py 123 data/fastqs/*.fq
```

random access to reads from fastqs

```sh
./benchmark_fastq_random_access.sh 3 data/fastqs/*.fq > benchmark_result_fastq_random_access.tsv
python3 make_fastq_ggplot2_matrix.py benchmark_result_fastq_random_access.tsv > matrix_fastq_random_access.tsv
```

### Sequence iteration

Iterating over sequences from fastas

```sh
./benchmark_fasta_sequence_iterate.sh 3 data/fastas/*.fa > benchmark_result_fasta_sequence_iterate.tsv
python3 make_fasta_ggplot2_matrix.py benchmark_result_fasta_sequence_iterate.tsv > matrix_fasta_sequence_iterate.tsv
```

Iterating over reads from fastqs

```sh
./benchmark_fastq_sequence_iterate.sh 3 data/fastqs/*.fq > benchmark_result_fastq_sequence_iterate.tsv
python3 make_fastq_ggplot2_matrix.py benchmark_result_fastq_sequence_iterate.tsv > matrix_fastq_sequence_iterate.tsv
```

### Plot

Recommended you use RStudio to plot, place the matrix file to the same folder of Rscripts. Prior to plotting, you should install ggplot2 and patchwork.

Select the work directory where matrix and rscripts located in

```r
setwd(dir.choose())
```

Plot for building index

```r
source('Figure_for_build_index.R')
```

Plot for random access

```r
source('Figure_for_random_access.R')
```

Plot for sequence iteration

```r
source('Figure_for_sequence_iteration.R')
```
