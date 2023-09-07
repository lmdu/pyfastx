Command line interface
======================

New in ``pyfastx`` 0.5.0

.. code:: bash

    $ pyfastx -h

    usage: pyfastx COMMAND [OPTIONS]

    A command line tool for FASTA/Q file manipulation

    optional arguments:
      -h, --help     show this help message and exit
      -v, --version  show program version number and exit

    Commands:

        index        build index for fasta/q file
        stat         show detailed statistics information of fasta/q file
        split        split fasta/q file into multiple files
        fq2fa        convert fastq file to fasta file
        subseq       get subsequences from fasta file by region
        sample       randomly sample sequences from fasta or fastq file
        extract      extract full sequences or reads from fasta/q file

Build index
-----------

New in ``pyfastx`` 0.6.10

.. code:: bash

    $ pyfastx index -h

    usage: pyfastx index [-h] [-f] fastx [fastx ...]

    positional arguments:
      fastx       fasta or fastq file, gzip support

    optional arguments:
      -h, --help  show this help message and exit
      -f, --full  build full index, base composition will be calculated

The --full option was used to count bases in FASTA/Q file and speedup calculation of GC content.

Show statistics information
---------------------------

.. code:: bash

    $ pyfastx stat -h

    usage: pyfastx stat [-h] fastx [fastx ...]

    positional arguments:
      fastx       input fasta or fastq file, gzip support

    optional arguments:
      -h, --help  show this help message and exit

For example:

.. code:: bash

    $ pyfastx info tests/data/*.fa*

    fileName        seqType seqCounts       totalBases         GC%   avgLen medianLen       maxLen  minLen  N50     L50
    protein.fa      protein        17             2265           -   133.24      80.0          419      23  263       4
    rna.fa              RNA         2              720      65.283    360.0     360.0          360     360  360       1
    test.fa             DNA       211            86262      43.529   408.82     386.0          821     118  516      66
    test.fa.gz          DNA       211            86262      43.529   408.82     386.0          821     118  516      66

seqType: sequence type (DNA, RNA, or protein); seqCounts: total sequence counts; totalBases: total number of bases; GC%: GC content; avgLen: average sequence length; medianLen: median sequence length; maxLen: maximum sequence length; minLen: minimum sequence length; N50: N50 length; L50: L50 sequence counts.

.. code:: bash

    $ pyfastx info tests/data/*.fq*

    fileName    readCounts  totalBases     GC%  avgLen  maxLen  minLen  maxQual minQual                     qualEncodingSystem
    test.fq            800      120000  66.175   150.0     150     150       70      35 Sanger Phred+33,Illumina 1.8+ Phred+33
    test.fq.gz         800      120000  66.175   150.0     150     150       70      35 Sanger Phred+33,Illumina 1.8+ Phred+33

readCounts: total read counts; totalBases: total number of bases; GC%: GC content; avgLen: average sequence length; maxLen: maximum sequence length; minLen: minimum sequence length; maxQual: maximum quality score; minQual: minimum quality score; qualEncodingSystem: quality encoding system.

Split FASTA/Q file
------------------

.. code:: bash

    $ pyfastx split -h

    usage: pyfastx split [-h] (-n int | -c int) [-o str] fastx

    positional arguments:
      fastx                 fasta or fastq file, gzip support

    optional arguments:
      -h, --help            show this help message and exit
      -n int                split a fasta/q file into N new files with even size
      -c int                split a fasta/q file into multiple files containing the same sequence counts
      -o str, --out-dir str
                            output directory, default is current folder

Convert FASTQ to FASTA file
---------------------------

.. code:: bash

    $ pyfastx fq2fa -h

    usage: pyfastx fq2fa [-h] [-o str] fastx

    positional arguments:
      fastx                 input fastq file, gzip support

    optional arguments:
      -h, --help            show this help message and exit
      -o str, --out-file str
                            output file, default: output to stdout

Get subsequence with region
---------------------------

.. code:: bash

    $ pyfastx subseq -h

    usage: pyfastx subseq [-h] [-r str | -b str] [-o str]
                          fastx [region [region ...]]

    positional arguments:
      fastx                 input fasta file, gzip support
      region                format is chr:start-end, start and end position is
                            1-based, multiple regions were separated by space

    optional arguments:
      -h, --help            show this help message and exit
      -r str, --region-file str
                            tab-delimited file, one region per line, both start
                            and end position are 1-based
      -b str, --bed-file str
                            tab-delimited BED file, 0-based start position and
                            1-based end position
      -o str, --out-file str
                            output file, default: output to stdout

Sample sequences
----------------

.. code:: bash

    $ pyfastx sample -h

    usage: pyfastx sample [-h] (-n int | -p float) [-s int] [--sequential-read]
                          [-o str]
                          fastx

    positional arguments:
      fastx                 fasta or fastq file, gzip support

    optional arguments:
      -h, --help            show this help message and exit
      -n int                number of sequences to be sampled
      -p float              proportion of sequences to be sampled, 0~1
      -s int, --seed int    random seed, default is the current system time
      --sequential-read     start sequential reading, particularly suitable for
                            sampling large numbers of sequences
      -o str, --out-file str
                            output file, default: output to stdout

Extract sequences
-----------------

New in ``pyfastx`` 0.6.10

.. code:: bash

    $ pyfastx extract -h

    usage: pyfastx extract [-h] [-l str] [--reverse-complement] [--out-fasta]
                           [-o str] [--sequential-read]
                           fastx [name [name ...]]

    positional arguments:
      fastx                 fasta or fastq file, gzip support
      name                  sequence name or read name, multiple names were
                            separated by space

    optional arguments:
      -h, --help            show this help message and exit
      -l str, --list-file str
                            a file containing sequence or read names, one name per
                            line
      --reverse-complement  output reverse complement sequence
      --out-fasta           output fasta format when extract reads from fastq,
                            default output fastq format
      -o str, --out-file str
                            output file, default: output to stdout
      --sequential-read     start sequential reading, particularly suitable for
                            extracting large numbers of sequences