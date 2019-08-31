Changelog
=========

Version 0.2.11 (2019-08-31)

- Support for Python 3.4

Version 0.2.10 (2019-08-28)
---------------------------

- Changed fseek and fread into gzseek and gzread
- Fixed sequence cache name comparision
- Fixed last sequence read error without line end
- Fixed subsequence slice error in normal FASTA file

Version 0.2.9 (2019-08-27)
--------------------------

- Fixed bad line calculate error
- Changed rewind to fseek for subsequence extraction

Version 0.2.8 (2019-08-26)
--------------------------

- Changed kseq.h library from li to attractivechaos
- Improved fasta parser

Version 0.2.7 (2019-08-26)
--------------------------

- Fixed no gzip index wrote to sqlite index file

Version 0.2.6 (2019-08-26)
--------------------------

- Optimized speed of gzip random access

Version 0.2.5 (2019-08-25)
--------------------------

- Fixed segmentation fault raised when loading gzip index
- Changed fasta object method get_seq to fetch

Version 0.2.4 (2019-08-25)
--------------------------

- Fixed fasta iter error after building new index

Version 0.2.3 (2019-08-24)
--------------------------

- Fixed fasta iter error when end of file is not \n

Version 0.2.2 (2019-07-19)
--------------------------

- Fixed identifier contain error

Version 0.2.1 (2019-07-15)
--------------------------

- Fixed sequence name always end with 0
- Fixed fasta iterable for flat fasta

Version 0.2.0 (2019-07-09)
--------------------------

- First release to PyPI
