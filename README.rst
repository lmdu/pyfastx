pyfastx
#######

.. image:: https://travis-ci.org/lmdu/pyfastx.svg?branch=master
   :target: https://travis-ci.org/lmdu/pyfastx
   :alt: Travis CI

.. image:: https://ci.appveyor.com/api/projects/status/7qeurb8wcl0bw993?svg=true
   :target: https://ci.appveyor.com/project/lmdu/pyfastx
   :alt: Appveyor CI

.. image:: https://readthedocs.org/projects/pyfastx/badge/?version=latest
   :target: https://pyfastx.readthedocs.io/en/latest/?badge=latest
   :alt: Readthedocs

.. image:: https://codecov.io/gh/lmdu/pyfastx/branch/master/graph/badge.svg
   :target: https://codecov.io/gh/lmdu/pyfastx
   :alt: Codecov

.. image:: https://coveralls.io/repos/github/lmdu/pyfastx/badge.svg?branch=master
   :target: https://coveralls.io/github/lmdu/pyfastx?branch=master
   :alt: Coveralls

.. image:: https://scan.coverity.com/projects/19696/badge.svg
   :target: https://scan.coverity.com/projects/lmdu-pyfastx
   :alt: Coverity

.. image:: https://img.shields.io/pypi/v/pyfastx.svg
   :target: https://pypi.org/project/pyfastx
   :alt: PyPI

.. image:: https://img.shields.io/pypi/wheel/pyfastx.svg
   :target: https://pypi.org/project/pyfastx
   :alt: Wheel

.. image:: https://img.shields.io/pypi/implementation/pyfastx
   :target: https://pypi.org/project/pyfastx
   :alt: PyPI - Implementation

.. image:: https://img.shields.io/pypi/pyversions/pyfastx.svg
   :target: https://pypi.org/project/pyfastx
   :alt: Pyver

.. image:: https://api.codacy.com/project/badge/Grade/80790fa30f444d9d9ece43689d512dae
   :target: https://www.codacy.com/manual/lmdu/pyfastx?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=lmdu/pyfastx&amp;utm_campaign=Badge_Grade
   :alt: Codacy

.. image:: https://img.shields.io/pypi/dm/pyfastx
   :target: https://pypi.org/project/pyfastx
   :alt: Downloads

.. image:: https://img.shields.io/pypi/l/pyfastx
   :target: https://pypi.org/project/pyfastx
   :alt: License

.. image:: https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat
   :target: http://bioconda.github.io/recipes/pyfastx/README.html
   :alt: Bioconda

*a robust python module for fast random access to sequences from plain and gzipped FASTA/Q file*

.. contents:: Table of Contents

Introduction
============

The ``pyfastx`` is a lightweight Python C extension that enables users to randomly access to sequences from plain and **gzipped** FASTA/Q files. This module aims to provide simple APIs for users to extract seqeunce from FASTA and reads from FASTQ by identifier and index number. The ``pyfastx`` will build indexes stored in a sqlite3 database file for random access to avoid consuming excessive amount of memory. In addition, the ``pyfastx`` can parse standard (*sequence is spread into multiple lines with same length*) and nonstandard (*sequence is spread into one or more lines with different length*) FASTA format. This module used `kseq.h <https://github.com/attractivechaos/klib/blob/master/kseq.h>`_ written by `@attractivechaos <https://github.com/attractivechaos>`_ in `klib <https://github.com/attractivechaos/klib>`_ project to parse plain FASTA/Q file and zran.c written by `@pauldmccarthy <https://github.com/pauldmccarthy>`_ in project `indexed_gzip <https://github.com/pauldmccarthy/indexed_gzip>`_ to index gzipped file for random access.

This project was heavily inspired by `@mdshw5 <https://github.com/mdshw5>`_'s project `pyfaidx <https://github.com/mdshw5/pyfaidx>`_ and `@brentp <https://github.com/brentp>`_'s project `pyfasta <https://github.com/brentp/pyfasta>`_.

Features
========

- Single file for the Python extension
- Lightweight, memory efficient for parsing FASTA/Q file
- Fast random access to sequences from ``gzipped`` FASTA/Q file
- Read sequences from FASTA file line by line
- Calculate N50 and L50 of sequences in FASTA file
- Calculate GC content and nucleotides composition
- Extract reverse, complement and antisense sequences
- Excellent compatibility, support for parsing nonstandard FASTA file
- Support for FASTQ quality score conversion
- Provide command line interface for splitting FASTA/Q file

Installation
============

Currently, ``pyfastx`` supports Python 3.5, 3.6, 3.7, 3.8. Make sure you have installed both `pip <https://pip.pypa.io/en/stable/installing/>`_ and Python before starting.

You can install ``pyfastx`` via the Python Package Index (PyPI)

::

    pip install pyfastx

Update ``pyfastx`` module

::

	pip install -U pyfastx

FASTA
=====

Read FASTA file
---------------

Read plain or gzipped FASTA file and build index, support for random access to FASTA.

.. code:: python

    >>> import pyfastx
    >>> fa = pyfastx.Fasta('test/data/test.fa.gz')
    >>> fa
    <Fasta> test/data/test.fa.gz contains 211 seqs

.. note::
    Building index may take some times. The time required to build index depends on the size of FASTA file. If index built, you can randomly access to any sequences in FASTA file. The index file can be reused to save time when you read seqeunces from FASTA file next time.

FASTA records iteration
-----------------------

The fastest way to iterate plain or gzipped FASTA file without building index, the iteration will return a tuple contains name and sequence.

.. code:: python

    >>> import pyfastx
    >>> for name, seq in pyfastx.Fasta('test/data/test.fa.gz', build_index=False):
    >>>     print(name, seq)

You can also iterate sequence object from FASTA object like this:

.. code:: python

    >>> import pyfastx
    >>> for seq in pyfastx.Fasta('test/data/test.fa.gz'):
    >>>     print(seq.name)
    >>>     print(seq.seq)
    >>>     print(seq.description)

Iteration with ``build_index=True`` (default) return sequence object which allows you to access attributions of sequence. New in pyfastx 0.6.3.


Get FASTA information
---------------------

.. code:: python

    >>> # get sequence counts in FASTA
    >>> len(fa)
    211

    >>> # get total sequence length of FASTA
    >>> fa.size
    86262

    >>> # get GC content of DNA sequence of FASTA
    >>> fa.gc_content
    43.529014587402344

    >>> # get GC skew of DNA sequences in FASTA
    >>> # New in pyfastx 0.3.8
    >>> fa.gc_skews
    0.004287730902433395

    >>> # get composition of nucleotides in FASTA
    >>> fa.composition
    {'A': 24534, 'C': 18694, 'G': 18855, 'T': 24179}

    >>> # get fasta type (DNA, RNA, or protein)
    >>> fa.type
    'DNA'

    >>> # check fasta file is gzip compressed
    >>> fa.is_gzip
    True

Get longest and shortest sequence
---------------------------------

New in ``pyfastx`` 0.3.0

.. code:: python

    >>> # get longest sequence
    >>> s = fa.longest
    >>> s
    <Sequence> JZ822609.1 with length of 821

    >>> s.name
    'JZ822609.1'

    >>> len(s)
    821

    >>> # get shortest sequence
    >>> s = fa.shortest
    >>> s
    <Sequence> JZ822617.1 with length of 118

    >>> s.name
    'JZ822617.1'

    >>> len(s)
    118

Calculate N50 and L50
---------------------

New in ``pyfastx`` 0.3.0

Calculate assembly N50 and L50, return (N50, L50), learn more about `N50,L50 <https://www.molecularecologist.com/2017/03/whats-n50/>`_

.. code:: python

	>>> # get FASTA N50 and L50
	>>> fa.nl(50)
	(516, 66)

	>>> # get FASTA N90 and L90
	>>> fa.nl(90)
	(231, 161)

	>>> # get FASTA N75 and L75
	>>> fa.nl(75)
	(365, 117)

Get sequence mean and median length
-----------------------------------

New in ``pyfastx`` 0.3.0

.. code:: python

	>>> # get sequence average length
	>>> fa.mean
	408

	>>> # get seqeunce median length
	>>> fa.median
	430

Get sequence counts
-------------------

New in ``pyfastx`` 0.3.0

Get counts of sequences whose length >= specified length

.. code:: python

	>>> # get counts of sequences with length >= 200 bp
	>>> fa.count(200)
	173

	>>> # get counts of sequences with length >= 500 bp
	>>> fa.count(500)
	70

Get subsequences
----------------

Subsequences can be retrieved from FASTA file by using a list of [start, end] coordinates

.. code:: python

    >>> # get subsequence with start and end position
    >>> interval = (1, 10)
    >>> fa.fetch('JZ822577.1', interval)
    'CTCTAGAGAT'

    >>> # get subsequences with a list of start and end position
    >>> intervals = [(1, 10), (50, 60)]
    >>> fa.fetch('JZ822577.1', intervals)
    'CTCTAGAGATTTTAGTTTGAC'

    >>> # get subsequences with reverse strand
    >>> fa.fetch('JZ822577.1', (1, 10), strand='-')
    'ATCTCTAGAG'

Key function
------------

New in ``pyfastx`` 0.5.1

Sometimes your fasta will have a long header which contains multiple identifiers and description, for example, ">JZ822577.1 contig1 cDNA library of flower petals in tree peony by suppression subtractive hybridization Paeonia suffruticosa cDNA, mRNA sequence". In this case, both "JZ822577.1" and "contig1" can be used as identifer. you can specify the key function to select one as identifier.

.. code:: python

	>>> #default use JZ822577.1 as identifier
	>>> #specify key_func to select contig1 as identifer
	>>> fa = pyfastx.Fasta('tests/data/test.fa.gz', key_func=lambda x: x.split()[1])
	>>> fa
	<Fasta> tests/data/test.fa.gz contains 211 seqs

Sequence
========

Get a sequence from FASTA
-------------------------

.. code:: python

    >>> # get sequence like a dictionary by identifier
    >>> s1 = fa['JZ822577.1']
    >>> s1
    <Sequence> JZ822577.1 with length of 333

    >>> # get sequence like a list by index
    >>> s2 = fa[2]
    >>> s2
    <Sequence> JZ822579.1 with length of 176

    >>> # get last sequence
    >>> s3 = fa[-1]
    >>> s3
    <Sequence> JZ840318.1 with length of 134

    >>> # check a sequence name weather in FASTA file
    >>> 'JZ822577.1' in fa
    True

Get sequence information
------------------------

.. code:: python

    >>> s = fa[-1]
    >>> s
    <Sequence> JZ840318.1 with length of 134

    >>> # get sequence order number in FASTA file
    >>> # New in pyfastx 0.3.7
    >>> s.id
    211

    >>> # get sequence name
    >>> s.name
    'JZ840318.1'

    >>> # get sequence description
    >>> # New in pyfastx 0.3.1
    >>> s.description
    'R283 cDNA library of flower petals in tree peony by suppression subtractive hybridization Paeonia suffruticosa cDNA, mRNA sequence'

    >>> # get sequence string
    >>> s.seq
    'ACTGGAGGTTCTTCTTCCTGTGGAAAGTAACTTGTTTTGCCTTCACCTGCCTGTTCTTCACATCAACCTTGTTCCCACACAAAACAATGGGAATGTTCTCACACACCCTGCAGAGATCACGATGCCATGTTGGT'

    >>> # get sequence raw string, New in pyfastx 0.6.3
    >>> print(s.raw)
    >JZ840318.1 R283 cDNA library of flower petals in tree peony by suppression subtractive hybridization Paeonia suffruticosa cDNA, mRNA sequence
    ACTGGAGGTTCTTCTTCCTGTGGAAAGTAACTTGTTTTGCCTTCACCTGCCTGTTCTTCACATCAACCTT
    GTTCCCACACAAAACAATGGGAATGTTCTCACACACCCTGCAGAGATCACGATGCCATGTTGGT

    >>> # get sequence length
    >>> len(s)
    134

    >>> # get GC content if dna sequence
    >>> s.gc_content
    46.26865768432617

    >>> # get nucleotide composition if dna sequence
    >>> s.composition
    {'A': 31, 'C': 37, 'G': 25, 'T': 41, 'N': 0}

Sequence slice
--------------

Sequence object can be sliced like a python string

.. code:: python

    >>> # get a sub seq from sequence
    >>> s = fa[-1]
    >>> ss = s[10:30]
    >>> ss
    <Sequence> JZ840318.1 from 11 to 30

    >>> ss.name
    'JZ840318.1:11-30'

    >>> ss.seq
    'CTTCTTCCTGTGGAAAGTAA'

    >>> ss = s[-10:]
    >>> ss
    <Sequence> JZ840318.1 from 125 to 134

    >>> ss.name
    'JZ840318.1:125-134'

    >>> ss.seq
    'CCATGTTGGT'


.. note::

	Slicing start and end coordinates are 0-based. Currently, pyfastx does not support an optional third ``step`` or ``stride`` argument. For example ``ss[::-1]``

Reverse and complement sequence
-------------------------------

.. code:: python

    >>> # get sliced sequence
    >>> fa[0][10:20].seq
    'GTCAATTTCC'

    >>> # get reverse of sliced sequence
    >>> fa[0][10:20].reverse
    'CCTTTAACTG'

    >>> # get complement of sliced sequence
    >>> fa[0][10:20].complement
    'CAGTTAAAGG'

    >>> # get reversed complement sequence, corresponding to sequence in antisense strand
    >>> fa[0][10:20].antisense
    'GGAAATTGAC'

Read sequence line by line
--------------------------

New in ``pyfastx`` 0.3.0

The sequence object can be iterated line by line as they appear in FASTA file.

.. code:: python

	>>> for line in fa[0]:
	... 	print(line)
	...
	CTCTAGAGATTACTTCTTCACATTCCAGATCACTCAGGCTCTTTGTCATTTTAGTTTGACTAGGATATCG
	AGTATTCAAGCTCATCGCTTTTGGTAATCTTTGCGGTGCATGCCTTTGCATGCTGTATTGCTGCTTCATC
	ATCCCCTTTGACTTGTGTGGCGGTGGCAAGACATCCGAAGAGTTAAGCGATGCTTGTCTAGTCAATTTCC
	CCATGTACAGAATCATTGTTGTCAATTGGTTGTTTCCTTGATGGTGAAGGGGCTTCAATACATGAGTTCC
	AAACTAACATTTCTTGACTAACACTTGAGGAAGAAGGACAAGGGTCCCCATGT

.. note::

    Sliced sequence (e.g. fa[0][10:50]) cannot be read line by line

Search for subsequence
----------------------

New in ``pyfastx`` 0.3.6

Search for subsequence from given sequence and get one-based start position of the first occurrence

.. code:: python

    >>> # search subsequence in sense strand
    >>> fa[0].search('GCTTCAATACA')
    262

    >>> # check subsequence weather in sequence
    >>> 'GCTTCAATACA' in fa[0]
    True

    >>> # search subsequence in antisense strand
    >>> fa[0].search('CCTCAAGT', '-')
    301

FASTQ
=====

New in ``pyfastx`` 0.4.0

Read FASTQ file
---------------

Read plain or gzipped file and build index, support for random access to reads from FASTQ.

.. code:: python

    >>> import pyfastx
    >>> fq = pyfastx.Fastq('tests/data/test.fq.gz')
    >>> fq
    <Fastq> tests/data/test.fq.gz contains 100 reads

FASTQ records iteration
-----------------------

The fastest way to parse plain or gzipped FASTQ file without building index, the iteration will return a tuple contains read name, seq and quality.

.. code:: python

    >>> import pyfastx
    >>> for name,seq,qual in pyfastx.Fastq('tests/data/test.fq.gz', build_index=False):
    >>>     print(name)
    >>>     print(seq)
    >>>     print(qual)

You can also iterate read object from FASTQ object like this:

.. code:: python

    >>> import pyfastx
    >>> for read in pyfastx.Fastq('test/data/test.fq.gz'):
    >>>     print(read.name)
    >>>     print(read.seq)
    >>>     print(read.qual)
    >>>     print(read.quali)

Iteration with ``build_index=True`` (default) return read object which allows you to access attribution of read. New in pyfastx 0.6.3.


Get FASTQ information
---------------------

.. code:: python

    >>> # get read counts in FASTQ
    >>> len(fq)
    800

    >>> # get total bases
    >>> fq.size
    120000

    >>> # get GC content of FASTQ file
    >>> fq.gc_content
    66.17471313476562

    >>> # get composition of bases in FASTQ
    >>> fq.composition
    {'A': 20501, 'C': 39705, 'G': 39704, 'T': 20089, 'N': 1}

    >>> # get phred which affects the quality score conversion
    >>> fq.phred
    33

    >>> # Guess fastq quality encoding system
    >>> # New in pyfastx 0.4.1
    >>> fq.encoding_type
    ['Sanger Phred+33', 'Illumina 1.8+ Phred+33']

Read
====

Get read from FASTQ
-------------------

.. code:: python

    >>> #get read like a dict by read name
    >>> r1 = fq['A00129:183:H77K2DMXX:1:1101:4752:1047']
    >>> r1
    <Read> A00129:183:H77K2DMXX:1:1101:4752:1047 with length of 150

    >>> # get read like a list by index
    >>> r2 = fq[10]
    >>> r2
    <Read> A00129:183:H77K2DMXX:1:1101:18041:1078 with length of 150

    >>> # get the last read
    >>> r3 = fq[-1]
    >>> r3
    <Read> A00129:183:H77K2DMXX:1:1101:31575:4726 with length of 150

    >>> # check a read weather in FASTQ file
    >>> 'A00129:183:H77K2DMXX:1:1101:4752:1047' in fq
    True

Get read information
--------------------

.. code:: python

    >>> r = fq[-10]
    >>> r
    <Read> A00129:183:H77K2DMXX:1:1101:1750:4711 with length of 150

    >>> # get read order number in FASTQ file
    >>> r.id
    791

    >>> # get read name
    >>> r.name
    'A00129:183:H77K2DMXX:1:1101:1750:4711'

    >>> # get read full header line, New in pyfastx 0.6.3
    >>> r.description
    '@A00129:183:H77K2DMXX:1:1101:1750:4711 1:N:0:CAATGGAA+CGAGGCTG'

    >>> # get read length
    >>> len(r)
    150

    >>> # get read sequence
    >>> r.seq
    'CGAGGAAATCGACGTCACCGATCTGGAAGCCCTGCGCGCCCATCTCAACCAGAAATGGGGTGGCCAGCGCGGCAAGCTGACCCTGCTGCCGTTCCTGGTCCGCGCCATGGTCGTGGCGCTGCGCGACTTCCCGCAGTTGAACGCGCGCTA'

    >>> # get raw string of read, New in pyfastx 0.6.3
    >>> print(r.raw)
    @A00129:183:H77K2DMXX:1:1101:1750:4711 1:N:0:CAATGGAA+CGAGGCTG
    CGAGGAAATCGACGTCACCGATCTGGAAGCCCTGCGCGCCCATCTCAACCAGAAATGGGGTGGCCAGCGCGGCAAGCTGACCCTGCTGCCGTTCCTGGTCCGCGCCATGGTCGTGGCGCTGCGCGACTTCCCGCAGTTGAACGCGCGCTA
    +
    FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FF,FFFFFFFFFFFFFFFFFFFFFFFFFF,F:FFFFFFFFF:

    >>> # get read quality ascii string
    >>> r.qual
    'FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FF,FFFFFFFFFFFFFFFFFFFFFFFFFF,F:FFFFFFFFF:'

    >>> # get read quality integer value, ascii - 33 or 64
    >>> r.quali
    [37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 25, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 25, 37, 37, 11, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 11, 37, 25, 37, 37, 37, 37, 37, 37, 37, 37, 37, 25]

    >>> # get read length
    >>> len(r)
    150

Identifiers
===========

Get identifiers
---------------

Get all identifiers of sequence as a list-like object.

.. code:: python

    >>> ids = fa.keys()
    >>> ids
    <Identifier> contains 211 identifiers

    >>> # get count of sequence
    >>> len(ids)
    211

    >>> # get identifier by index
    >>> ids[0]
    'JZ822577.1'

    >>> # check identifier where in fasta
    >>> 'JZ822577.1' in ids
    True

    >>> # iter identifiers
    >>> for name in ids:
    >>>     print(name)

    >>> # convert to a list
    >>> list(ids)

Sort identifiers
----------------

Sort identifiers by sequence id, name, or length for iteration

New in ``pyfastx`` 0.5.0

.. code:: python

    >>> # sort identifiers by length with descending order
    >>> for name in ids.sort(key='length', reverse=True):
    >>>     print(name)

    >>> # sort identifiers by name with ascending order
    >>> for name in ids.sort(key='name'):
    >>>     print(name)

    >>> # sort identifiers by id with descending order
    >>> for name in ids.sort(key='id', reverse=True)
    >>>     print(name)

Filter identifiers
------------------

Filter identifiers by sequence length and name

New in ``pyfastx`` 0.5.10

.. code:: python

    >>> # get identifiers with length > 600
    >>> ids.filter(ids > 600)
    <Identifier> contains 48 identifiers

    >>> # get identifiers with length >= 500 and <= 700
    >>> ids.filter(ids>=500, ids<=700)
    <Identifier> contains 48 identifiers

    >>> # get identifiers with length > 500 and < 600
    >>> ids.filter(500<ids<600)
    <Identifier> contains 22 identifiers

    >>> # get identifiers contain JZ8226
    >>> ids.filter(ids % 'JZ8226')
    <Identifier> contains 90 identifiers

    >>> # get identifiers contain JZ8226 with length > 550
    >>> ids.filter(ids % 'JZ8226', ids>550)
    <Identifier> contains 17 identifiers

    >>> # clear sort order and filters
    >>> ids.reset()
    <Identifier> contains 211 identifiers

    >>> # list a filtered result
    >>> ids.filter(ids % 'JZ8226', ids>730)
    >>> list(ids)
    ['JZ822609.1', 'JZ822650.1', 'JZ822664.1', 'JZ822699.1']

    >>> # list a filtered result with sort order
    >>> ids.filter(ids % 'JZ8226', ids>730).sort('length', reverse=True)
    >>> list(ids)
    ['JZ822609.1', 'JZ822699.1', 'JZ822664.1', 'JZ822650.1']

    >>> ids.filter(ids % 'JZ8226', ids>730).sort('name', reverse=True)
    >>> list(ids)
    ['JZ822699.1', 'JZ822664.1', 'JZ822650.1', 'JZ822609.1']

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

Drawbacks
=========

If you intensively check sequence names exists in FASTA file using ``in`` operator on FASTA object like:

.. code:: python

	>>> fa = pyfastx.Fasta('tests/data/test.fa.gz')
	>>> # Suppose seqnames has 100000 names
	>>> for seqname in seqnames:
	>>>     if seqname in fa:
	>>>	        do something

This will take a long time to finish. Becuase, pyfastx does not load the index into memory, the ``in`` operating is corresponding to sql query existence from index database. The faster alternative way to do this is:

.. code:: python

	>>> fa = pyfastx.Fasta('tests/data/test.fa.gz')
	>>> # load all sequence names into a set object
	>>> all_names = set(fa.keys())
	>>> for seqname in seqnames:
	>>>     if seqname in all_names:
	>>>	        do something

Testing
=======

The ``pyfaidx`` module was used to test ``pyfastx``. To run the tests:

::

	$ python setup.py test

Acknowledgements
================

`kseq.h <https://github.com/attractivechaos/klib/blob/master/kseq.h>`_ and `zlib <https://www.zlib.net/>`_ was used to parse FASTA format. `Sqlite3 <https://www.sqlite.org/index.html>`_ was used to store built indexes. ``pyfastx`` can randomly access to sequences from gzipped FASTA file mainly attributed to `indexed_gzip <https://github.com/pauldmccarthy/indexed_gzip>`_.
