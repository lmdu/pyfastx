Command line interface
======================

New in ``pyfastx`` 0.5.0

.. code:: bash

    $ pyfastx -h
    
    usage: pyfastx COMMAND [OPTIONS]

    A command line tool for FASTA/Q file manipulation

    optional arguments:
      -h, --help     show this help message and exit
      -v, --version  show program's version number and exit

    Commands:

        index        build index for FASTA or FASTQ file
        info         show detailed statistics information of FASTA/Q file
        split        split fasta file into multiple files
        fq2fa        convert fastq file to fasta file
        subseq       get subseqence from fasta file by id or name with region
        sample       randomly sample sequences from fasta or fastq file
        extract      extract sequences or reads from fasta or fastq file

Build index
-----------

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

    $ pyfastx info -h

    usage: pyfastx info [-h] fastx [fastx ...]

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

..code:: bash

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
      -n int                split a fa/q file into N new files with even size
      -c int                split a fa/q file into multiple files with the same
                            sequence counts
      -o str, --outdir str  output directory, default is current folder

Convert FASTQ to FASTA file
---------------------------

.. code:: bash

    $ pyfastx fq2fa -h

    usage: pyfastx fq2fa [-h] [-o str] fastx

    positional arguments:
      fastx                 input fastq file, gzip support

    optional arguments:
      -h, --help            show this help message and exit
      -o str, --outfile str
                            output file, default: output to stdout

Get subsequence with region
---------------------------

.. code:: bash

    $ pyfastx subseq -h

    usage: pyfastx subseq [-h] (--id int | --chr str) [-r str] fastx

    positional arguments:
      fastx                 input fasta file, gzip support

    optional arguments:
      -h, --help            show this help message and exit
      --id int              sequence id number in fasta file
      --chr str             sequence name
      -r str, --region str  one-based slice region, e.g. 10:20

Sample sequences
----------------

.. code:: bash

    $ pyfastx sample -h

    usage: pyfastx sample [-h] (-n int | -p float) [-o str] fastx

    positional arguments:
      fastx                 fasta or fastq file, gzip support

    optional arguments:
      -h, --help            show this help message and exit
      -n int                number of sequences to be sampled
      -p float              proportion of sequences to be sampled, 0~1
      -o str, --outfile str
                            output file, default: output to stdout

Extract sequences
-----------------

.. code:: bash

    $ pyfastx extract -h

    usage: pyfastx extract [-h] (--ids int or str | --names str) [--outfas]
                           [-o str]
                           fastx

    positional arguments:
      fastx                 fasta or fastq file, gzip support

    optional arguments:
      -h, --help            show this help message and exit
      --ids int or str      extract sequences by id number, the value can be one
                            integer to get one sequence, a range (e.g. 5-10) or a
                            comma seperated list (e.g. 3,5,8) to get multiple
                            sequences
      --names str           extract sequences by name, the value can be one name
                            to get one sequence, a comma seperated list (e.g.
                            seq1,seq5,seq9) or a file contains names (one name per
                            line) to get multiple sequences
      --outfas              output fasta format when input file is fastq format,
                            default output fastq format
      -o str, --outfile str
                            output file, default: output to stdout