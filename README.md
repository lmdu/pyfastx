# pyfastx

*a robust python module for fast random access to FASTA sequences*

* [About](#about)
* [Installation](#installation)
* [Usage](#usage)
* [Acknowledgements](#acknowledgements)
* [License](#license)



## About

The pyfastx is a lightweight Python C extension that enables you to randomly access FASTA sequences  in flat text file, even in gzip compressed file. This module uses [kseq.h](http://lh3lh3.users.sourceforge.net/kseq.shtml) written by Heng Li to parse FASTA file and zran.c written by Paul McCarthy in [indexed_gzip](https://github.com/pauldmccarthy/indexed_gzip) project to index gzipped file.



## Installation

`pyfastx`  is available on [PyPi](https://pypi.org), to install, simply type:

```sh
pip install pyfastx
```



## Usage

You can use the pyfastx module directly:

```python
>>> import pyfastx
>>> fa = pyfastx.Fasta('test.fa.gz')
>>> fa
Fasta('test.fa.gz')
```

