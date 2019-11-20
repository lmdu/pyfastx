Changelog
=========

Version 0.5.10 (2019-11-20)
---------------------------

- Added identifier filter function
- Remove tp_new for Read, Sequence and Identifier
- Fixed module method error

Version 0.5.9 (2019-11-17)
--------------------------

- Added get longest and shortest sequence object
- Added composition argument to speedup getting GC content
- Added memory index to keep index in memory rather than local file
- Fixed command line error
- Changed sqlite to higher version
- Removed journal_mode OFF
- Speedup index building

Version 0.5.8 (2019-11-10)
--------------------------

- Fixed fasta NL function parameter check
- Fixed read id error when fastq iteration

Version 0.5.7 (2019-11-09)
--------------------------

- Fixed SystemError caused caused by Python 2.7 seperated int and long type
- Fixed String type check on Python 2.7
- Fixed objects memory deallocation

Version 0.5.6 (2019-11-08)
--------------------------

- Optimized random access from plain file
- Reduced memory consumption

Version 0.5.5 (2019-11-07)
--------------------------

- Added Support for IUPAC code complement
- Speedup reverse complement
- Speedup space removing and uppercase


Version 0.5.4 (2019-11-04)
--------------------------

- Added guess fasta type (DNA, RNA, protein)
- Added support for calculating protein sequence composition
- Optimized the speed of index building
- Calculate sequence composition when get gc content or composition
- Fixed char return in python 2.7

Version 0.5.3 (2019-10-23)
--------------------------

- Added support for coverting fastq to fasta
- Updated command line interface docs
- Fixed command line entry points

Version 0.5.2 (2019-10-18)
--------------------------

- Fixed command line interface running error

Version 0.5.1 (2019-10-17)
--------------------------

- Added key function for custom sequence identifier
- Optimized speed of fasta indexing
- Fixed bool args parsing error in py2.7

Version 0.5.0 (2019-10-13)
--------------------------

- Added support for python 2.7 and 3.4
- Added command line tool to manipulate fasta and fastq file
- Added gzip attribute to fasta and fastq object to check whether compressed
- Added sort function for identifier object
- Fixed python bool argument parsing error caused by uint16_t
- Fixed identifier sort key initialization

Version 0.4.1 (2019-10-05)
--------------------------

- Fixed fastq quality encoding system guesser
- Fixed gzip index insertion error

Version 0.4.0 (2019-09-29)
--------------------------

- Added support for parsing FASTQ
- Added random access to reads from FASTQ

Version 0.3.10 (2019-09-27)
---------------------------

- Fixed GC skew exception caused by mixing unsigned with signed for division  

Version 0.3.9 (2019-09-26)
--------------------------

- Fixed sequence read line by line error
- Fixed last sequence build index error when fasta file ended without \n
- Fixed GC skew error

Version 0.3.8 (2019-09-25)
--------------------------

- Fixed large offset became negative error
- Fixed slice step
- Fixed uncorrect median length
- Fixed strand compare error
- Added GC skew calculation
- Updated test script

Version 0.3.7 (2019-09-24)
--------------------------

- Changed int type to standard type
- Added support for processing large fasta file
- Added id number for each sequence
- Fixed SQL fetch error
- Used 50 as default value of nl to calculate N50 and L50

Version 0.3.6 (2019-09-20)
--------------------------

- Added support for searching subsequence from a sequence
- Added support for checking subsequence weather in a sequence
- Fixed gzip index import error
- Fixed subsequence parent length for full sequence extraction

Version 0.3.5 (2019-09-08)
--------------------------

- Fixed unicode error caused by sqlite3_finalize 

Version 0.3.4 (2019-09-07)
--------------------------

- Fixed seq description unicode string error

Version 0.3.3 (2019-09-07)
--------------------------

- Fixed sequence description encoding error
 
Version 0.3.2 (2019-09-07)
--------------------------

Deleted

Version 0.3.1 (2019-09-07)
--------------------------

- Added support for geting sequence description

Version 0.3.0 (2019-09-07)
--------------------------

- Added read sequence from fasta file line by line
- Added support for calculating assembly N50 and L50
- Added support for calculating median and average length
- Added support for getting longest and shortest sequence
- Added support for calculating counts of sequence
- removed support for Python34

Version 0.2.11 (2019-08-31)
---------------------------

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
