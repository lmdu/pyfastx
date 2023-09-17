pyfastx
#######

.. image:: https://github.com/lmdu/pyfastx/actions/workflows/main.yml/badge.svg
   :target: https://github.com/lmdu/pyfastx/actions/workflows/main.yml
   :alt: Action

.. image:: https://readthedocs.org/projects/pyfastx/badge/?version=latest
   :target: https://pyfastx.readthedocs.io/en/latest/?badge=latest
   :alt: Readthedocs

.. image:: https://codecov.io/gh/lmdu/pyfastx/branch/master/graph/badge.svg
   :target: https://codecov.io/gh/lmdu/pyfastx
   :alt: Codecov

.. image:: https://img.shields.io/pypi/v/pyfastx.svg
   :target: https://pypi.org/project/pyfastx
   :alt: PyPI

.. image:: https://img.shields.io/pypi/implementation/pyfastx
   :target: https://pypi.org/project/pyfastx
   :alt: Language

.. image:: https://img.shields.io/pypi/pyversions/pyfastx.svg
   :target: https://pypi.org/project/pyfastx
   :alt: Pyver

.. image:: https://img.shields.io/pypi/wheel/pyfastx.svg
   :target: https://pypi.org/project/pyfastx
   :alt: Wheel

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

**Citation:** 
`Lianming Du, Qin Liu, Zhenxin Fan, Jie Tang, Xiuyue Zhang, Megan Price, Bisong Yue, Kelei Zhao. Pyfastx: a robust Python package for fast random access to sequences from plain and gzipped FASTA/Q files. Briefings in Bioinformatics, 2021, 22(4):bbaa368 <https://doi.org/10.1093/bib/bbaa368>`_.

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

Currently, ``pyfastx`` supports Python 3.6, 3.7, 3.8, 3.9, 3.10, 3.11. Make sure you have installed both `pip <https://pip.pypa.io/en/stable/installing/>`_ and Python before starting.

You can install ``pyfastx`` via the Python Package Index (PyPI)

::

    pip install pyfastx

Update ``pyfastx`` module

::

	pip install -U pyfastx

FASTX
=====

New in ``pyfastx`` 0.8.0.

Pyfastx provide a simple and fast python binding for kseq.h to iterate over sequences or reads in fasta/q file. The FASTX object will automatically detect the input sequence format (fasta or fastq) to return different tuple.

FASTA sequences iteration
-------------------------

When iterating over sequences on FASTX object, a tuple ``(name, seq)`` will be returned.

.. code:: python

    >>> fa = pyfastx.Fastx('tests/data/test.fa.gz')
    >>> for name,seq in fa:
    >>>     print(name)
    >>>     print(seq)

    >>> #always output uppercase sequence
    >>> for item in pyfastx.Fastx('tests/data/test.fa', uppercase=True):
    >>>     print(item)

    >>> #Manually specify sequence format
    >>> for item in pyfastx.Fastx('tests/data/test.fa', format="fasta"):
    >>>     print(item)

If you want the sequence comment, you can set comment to True, New in ``pyfastx`` 0.9.0.

.. code:: python

    >>> fa = pyfastx.Fastx('tests/data/test.fa.gz', comment=True)
    >>> for name,seq,comment in fa:
    >>>     print(name)
    >>>     print(seq)
    >>>     print(comment)

The comment is the content of header line after the first white space or tab character.

FASTQ reads iteration
---------------------

When iterating over reads on FASTX object, a tuple ``(name, seq, qual)`` will be returned.

.. code:: python

    >>> fq = pyfastx.Fastx('tests/data/test.fq.gz')
    >>> for name,seq,qual in fq:
    >>>     print(name)
    >>>     print(seq)
    >>>     print(qual)

If you want the read comment, you can set comment to True, New in ``pyfastx`` 0.9.0.

.. code:: python

    >>> fq = pyfastx.Fastx('tests/data/test.fq.gz', comment=True)
    >>> for name,seq,qual,comment in fq:
    >>>     print(name)
    >>>     print(seq)
    >>>     print(qual)
    >>>     print(comment)

The comment is the content of header line after the first white space or tab character.

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

FastaKeys
=========

New in ``pyfastx`` 0.8.0. We have changed ``Identifier`` object to ``FastaKeys`` object.

Get keys
--------------

Get all names of sequence as a list-like object.

.. code:: python

    >>> ids = fa.keys()
    >>> ids
    <FastaKeys> contains 211 keys

    >>> # get count of sequence
    >>> len(ids)
    211

    >>> # get key by index
    >>> ids[0]
    'JZ822577.1'

    >>> # check key whether in fasta
    >>> 'JZ822577.1' in ids
    True

    >>> # iterate over keys
    >>> for name in ids:
    >>>     print(name)

    >>> # convert to a list
    >>> list(ids)

Sort keys
----------------

Sort keys by sequence id, name, or length for iteration

New in ``pyfastx`` 0.5.0

.. code:: python

    >>> # sort keys by length with descending order
    >>> for name in ids.sort(by='length', reverse=True):
    >>>     print(name)

    >>> # sort keys by name with ascending order
    >>> for name in ids.sort(by='name'):
    >>>     print(name)

    >>> # sort keys by id with descending order
    >>> for name in ids.sort(by='id', reverse=True)
    >>>     print(name)

Filter keys
------------------

Filter keys by sequence length and name

New in ``pyfastx`` 0.5.10

.. code:: python

    >>> # get keys with length > 600
    >>> ids.filter(ids > 600)
    <FastaKeys> contains 48 keys

    >>> # get keys with length >= 500 and <= 700
    >>> ids.filter(ids>=500, ids<=700)
    <FastaKeys> contains 48 keys

    >>> # get keys with length > 500 and < 600
    >>> ids.filter(500<ids<600)
    <FastaKeys> contains 22 keys

    >>> # get keys contain JZ8226
    >>> ids.filter(ids % 'JZ8226')
    <FastaKeys> contains 90 keys

    >>> # get keys contain JZ8226 with length > 550
    >>> ids.filter(ids % 'JZ8226', ids>550)
    <FastaKeys> contains 17 keys

    >>> # clear sort order and filters
    >>> ids.reset()
    <FastaKeys> contains 211 keys

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

    >>> # New in pyfastx 0.6.10
    >>> # get average length of reads
    >>> fq.avglen
    150.0

    >>> # get maximum lenth of reads
    >>> fq.maxlen
    150

    >>> # get minimum length of reas
    >>> fq.minlen
    150

    >>> # get maximum quality score
    >>> fq.maxqual
    70

    >>> # get minimum quality score
    >>> fq.minqual
    35

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

FastqKeys
=========

New in ``pyfastx`` 0.8.0.

Get fastq keys
---------------

Get all names of read as a list-like object.

.. code:: python

    >>> ids = fq.keys()
    >>> ids
    <FastqKeys> contains 800 keys

    >>> # get count of read
    >>> len(ids)
    800

    >>> # get key by index
    >>> ids[0]
    'A00129:183:H77K2DMXX:1:1101:6804:1031'

    >>> # check key whether in fasta
    >>> 'A00129:183:H77K2DMXX:1:1101:14416:1031' in ids
    True

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

Show statistics information
---------------------------

.. code:: bash

    $ pyfastx stat -h

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
      fastx                 fastq file, gzip support

    optional arguments:
      -h, --help            show this help message and exit
      -o str, --out-file str
                            output file, default: output to stdout

Get subsequence with region
---------------------------

.. code:: bash

    $ pyfastx subseq -h

    usage: pyfastx subseq [-h] [-r str | -b str] [-o str] fastx [region [region ...]]

    positional arguments:
      fastx                 input fasta file, gzip support
      region                format is chr:start-end, start and end position is 1-based, multiple names were separated by space

    optional arguments:
      -h, --help            show this help message and exit
      -r str, --region-file str
                            tab-delimited file, one region per line, both start and end position are 1-based
      -b str, --bed-file str
                            tab-delimited BED file, 0-based start position and 1-based end position
      -o str, --out-file str
                            output file, default: output to stdout

Sample sequences
----------------

.. code:: bash

    $ pyfastx sample -h

    usage: pyfastx sample [-h] (-n int | -p float) [-s int] [--sequential-read] [-o str] fastx

    positional arguments:
      fastx                 fasta or fastq file, gzip support

    optional arguments:
      -h, --help            show this help message and exit
      -n int                number of sequences to be sampled
      -p float              proportion of sequences to be sampled, 0~1
      -s int, --seed int    random seed, default is the current system time
      --sequential-read     start sequential reading, particularly suitable for sampling large numbers of sequences
      -o str, --out-file str
                            output file, default: output to stdout

Extract sequences
-----------------

New in ``pyfastx`` 0.6.10

.. code:: bash

    $ pyfastx extract -h

    usage: pyfastx extract [-h] [-l str] [--reverse-complement] [--out-fasta] [-o str] [--sequential-read]
                           fastx [name [name ...]]

    positional arguments:
      fastx                 fasta or fastq file, gzip support
      name                  sequence name or read name, multiple names were separated by space

    optional arguments:
      -h, --help            show this help message and exit
      -l str, --list-file str
                            a file containing sequence or read names, one name per line
      --reverse-complement  output reverse complement sequence
      --out-fasta           output fasta format when extract reads from fastq, default output fastq format
      -o str, --out-file str
                            output file, default: output to stdout
      --sequential-read     start sequential reading, particularly suitable for extracting large numbers of sequences

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

The ``pyfaidx`` module was used to test ``pyfastx``. First, make sure you have a suitable version installed::

    pip install pyfastx

To test pyfastx, you should also install pyfaidx 0.5.8::

    pip install pyfaidx==0.5.8

Then, to run the tests::

	$ python setup.py test

Acknowledgements
================

`kseq.h <https://github.com/attractivechaos/klib/blob/master/kseq.h>`_ and `zlib <https://www.zlib.net/>`_ was used to parse FASTA format. `Sqlite3 <https://www.sqlite.org/index.html>`_ was used to store built indexes. ``pyfastx`` can randomly access to sequences from gzipped FASTA file mainly attributed to `indexed_gzip <https://github.com/pauldmccarthy/indexed_gzip>`_.
