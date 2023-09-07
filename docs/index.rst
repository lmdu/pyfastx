.. pyfastx documentation master file, created by
   sphinx-quickstart on Mon Jun 24 21:13:35 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to pyfastx's documentation!
===================================

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

.. image:: https://img.shields.io/pypi/pyversions/pyfastx.svg
   :target: https://pypi.org/project/pyfastx
   :alt: Pyver

.. image:: https://img.shields.io/pypi/wheel/pyfastx.svg
   :target: https://pypi.org/project/pyfastx
   :alt: Wheel

.. image:: https://api.codacy.com/project/badge/Grade/80790fa30f444d9d9ece43689d512dae
   :target: https://www.codacy.com/manual/lmdu/pyfastx?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=lmdu/pyfastx&amp;utm_campaign=Badge_Grade
   :alt: Codacy

.. image:: https://img.shields.io/pypi/implementation/pyfastx
   :target: https://pypi.org/project/pyfastx
   :alt: Language

.. image:: https://img.shields.io/pypi/dm/pyfastx
   :target: https://pypi.org/project/pyfastx
   :alt: Downloads

.. image:: https://img.shields.io/pypi/l/pyfastx
   :target: https://pypi.org/project/pyfastx
   :alt: License

.. image:: https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat
   :target: http://bioconda.github.io/recipes/pyfastx/README.html
   :alt: Bioconda

The ``pyfastx`` is a lightweight Python C extension that enables users to randomly access to sequences from plain and **gzipped** FASTA/Q files. This module aims to provide simple APIs for users to extract sequence from FASTA and reads from FASTQ by identifier and index number. The ``pyfastx`` will build indexes stored in a sqlite3 database file for random access to avoid consuming excessive amount of memory. In addition, the ``pyfastx`` can parse standard (*sequence is spread into multiple lines with same length*) and nonstandard (*sequence is spread into one or more lines with different length*) FASTA format. This module used `kseq.h <https://github.com/attractivechaos/klib/blob/master/kseq.h>`_ written by `@attractivechaos <https://github.com/attractivechaos>`_ in `klib <https://github.com/attractivechaos/klib>`_ project to parse plain FASTA/Q file and zran.c written by `@pauldmccarthy <https://github.com/pauldmccarthy>`_ in project `indexed_gzip <https://github.com/pauldmccarthy/indexed_gzip>`_ to index gzipped file for random access.

This project was heavily inspired by `@mdshw5 <https://github.com/mdshw5>`_'s project `pyfaidx <https://github.com/mdshw5/pyfaidx>`_ and `@brentp <https://github.com/brentp>`_'s project `pyfasta <https://github.com/brentp/pyfasta>`_.

**Features**

- Single file for the Python extension
- Lightweight, memory efficient for parsing FASTA file
- Fast random access to sequences from **gzipped** FASTA file
- Read sequences from FASTA file line by line
- Calculate assembly N50 and L50
- Calculate GC content and nucleotides composition
- Extract reverse, complement and antisense sequence
- Excellent compatibility, support for parsing nonstandard FASTA file
- Support for random access reads from FASTQ file

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   installation
   usage
   commandline
   advance
   drawbacks
   changelog
   api_reference
   acknowledgements


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
