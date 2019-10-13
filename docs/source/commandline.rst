Command line interface
======================

New in ``pyfastx`` 0.5.0

.. code:: bash

    $ python3 pyfastx.py -h
    
    usage: pyfastx COMMAND [OPTIONS]

    A command line tool for FASTA/Q file manipulation

    optional arguments:
      -h, --help     show this help message and exit
      -v, --version  show program's version number and exit

    Commands:

        info         Get detailed statistics information of FASTA/Q file
        split        Split fasta file into multiple files

Show statistics information
---------------------------

.. code:: bash

    $ python3 pyfastx.py info -h

    usage: pyfastx info [-h] fastx

    positional arguments:
      fastx       input fasta or fastq file, gzip compressed support

    optional arguments:
      -h, --help  show this help message and exit

Split FASTA/Q file
------------------

.. code:: bash

    $ python3 pyfastx.py split -h

    usage: pyfastx split [-h] (-n FILE_NUM | -c SEQ_COUNT) [-o OUT_DIR] [-g] fastx

    positional arguments:
      fastx                 input fasta or fastq file, gzip compressed support

    optional arguments:
      -h, --help            show this help message and exit
      -n FILE_NUM, --file_num FILE_NUM
                            split a fasta or fastq file into N new files with even
                            size
      -c SEQ_COUNT, --seq_count SEQ_COUNT
                            split a fasta or fastq file into multiple files with
                            the same sequence counts
      -o OUT_DIR, --out_dir OUT_DIR
                            output directory, default is current folder
      -g, --gzip_compress   use gzip to compress output files
