pyfastx
=======

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

.. image:: https://img.shields.io/pypi/v/pyfastx.svg
   :target: https://pypi.org/project/pyfastx
   :alt: PyPI Version

.. image:: https://img.shields.io/pypi/pyversions/pyfastx.svg
   :target: https://pypi.org/project/pyfastx
   :alt: Python Version

.. image:: https://img.shields.io/pypi/wheel/pyfastx.svg
   :target: https://pypi.org/project/pyfastx
   :alt: Wheel

*a robust python module for fast random access to sequences from plain and gzipped FASTA file*

.. contents:: Table of Contents
	:depth: 2

About
-----

The ``pyfastx`` is a lightweight Python C extension that enables users to randomly access to sequences from plain and **gzipped** FASTA files. This module aims to provide simple APIs for users to extract seqeunce from FASTA by identifier and index number. The ``pyfastx`` will build indexes stored in a sqlite3 database file for random access to avoid consuming excessive amount of memory. In addition, the ``pyfastx`` can parse standard (*sequence spread into multiple lines with same length*) and nonstandard (*lines with different length*) FASTA format. This module used `kseq.h <https://github.com/attractivechaos/klib/blob/master/kseq.h>`_ written by `@attractivechaos <https://github.com/attractivechaos>`_ in `klib <https://github.com/attractivechaos/klib>`_ project to parse plain FASTA file and zran.c written by `@pauldmccarthy <https://github.com/pauldmccarthy>`_ in project `indexed_gzip <https://github.com/pauldmccarthy/indexed_gzip>`_ to index gzipped file for random access.

This project was heavily inspired by `@mdshw5 <https://github.com/mdshw5>`_'s project `pyfaidx <https://github.com/mdshw5/pyfaidx>`_ and `@brentp <https://github.com/brentp>`_'s project `pyfasta <https://github.com/brentp/pyfasta>`_.

Installation
------------

Make sure you have both `pip <https://pip.pypa.io/en/stable/installing/>`_ and at least version 3.5 of Python before starting.

You can install ``pyfastx`` via the Python Package Index (PyPI)

::

    pip install pyfastx

Update ``pyfastx`` module

::

	pip install -U pyfastx

Usage
-----

Read FASTA file
^^^^^^^^^^^^^^^

The fastest way to parse flat or gzipped FASTA file without building index.

.. code:: python

    >>> import pyfastx
    >>> for name, seq in pyfastx.Fasta('test/data/test.fa.gz', build_index=False):
    >>>     print(name, seq)

Read flat or gzipped FASTA file and build index, support for random access to FASTA.

.. code:: python

    >>> import pyfastx
    >>> fa = pyfastx.Fasta('test/data/test.fa.gz')
    >>> fa
    <Fasta> test/data/test.fa.gz contains 211 seqs

.. note::

	Note: Building index may take some times. The time required to build index depends on the size of FASTA file. If index built, you can randomly access to any sequences in FASTA file.

Get FASTA information
^^^^^^^^^^^^^^^^^^^^^

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

    >>> # get composition of nucleotides in FASTA
    >>> fa.composition
    {'A': 24534, 'C': 18694, 'G': 18855, 'T': 24179, 'N': 0}

Get sequence from FASTA
^^^^^^^^^^^^^^^^^^^^^^^

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
^^^^^^^^^^^^^^^^^^^^^^^^

.. code:: python

    >>> s = fa[-1]
    >>> s
    <Sequence> JZ840318.1 with length of 134

    >>> # get sequence name
    >>> s.name
    'JZ840318.1'

    >>> # get sequence string
    >>> s.seq
    'ACTGGAGGTTCTTCTTCCTGTGGAAAGTAACTTGTTTTGCCTTCACCTGCCTGTTCTTCACATCAACCTTGTTCCCACACAAAACAATGGGAATGTTCTCACACACCCTGCAGAGATCACGATGCCATGTTGGT'

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
^^^^^^^^^^^^^^

Sequence object can be sliced like a python string

.. code:: python

    >>> # get a sub seq from sequence
    >>> ss = seq[10:30]
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

	Note: Slicing start and end coordinates are 0-based. Currently, pyfastx does not support an optional third ``step`` or ``stride`` argument. For example ``ss[::-1]``

Reverse and complement sequence
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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

Get subsequences
^^^^^^^^^^^^^^^^

Subseuqneces can be retrieved from FASTA file by using a list of [start, end] coordinates

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

Get identifiers
^^^^^^^^^^^^^^^

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

Testing
-------

The ``pyfaidx`` module was used to test ``pyfastx``. To run the tests:

::

	$ python setup.py test

Acknowledgements
----------------

`kseq.h <https://github.com/attractivechaos/klib/blob/master/kseq.h>`_ and `zlib <https://www.zlib.net/>`_ was used to parse FASTA format. `Sqlite3 <https://www.sqlite.org/index.html>`_ was used to store built indexes. ``pyfastx`` can randomly access to sequences from gzipped FASTA file mainly attributed to `indexed_gzip <https://github.com/pauldmccarthy/indexed_gzip>`_.
