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

        info         show detailed statistics information of FASTA/Q file
        split        split fasta file into multiple files
        fq2fa        Convert fastq file to fasta file
        subseq       Get subseqence from fasta file by id or name with region
        sample       randomly sample sequences from fasta or fastq file

Show statistics information
---------------------------

.. code:: bash

    $ pyfastx info -h

    usage: pyfastx info [-h] fastx

    positional arguments:
      fastx       input fasta or fastq file, gzip support

    optional arguments:
      -h, --help  show this help message and exit

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