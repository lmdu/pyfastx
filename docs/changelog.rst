Changelog
=========

Version 2.0.1 (2023-09-18)
--------------------------

- Speedup the gzip index writing to index file

Version 2.0.0 (2023-09-05)
--------------------------

- Added support for file name with wide char
- Added support for specifying index file path
- Added support for more characters in DNA sequence
- Added reverse complement function for DNA conversion
- Improved the performance of kseq library
- Optimized gzip index importing and saving without temp file
- Fixed segmentation fault when using sequence composition
- Fixed memory leak in Fastq read quality integer
- Fixed zlib download url broken error when building

Version 1.1.0 (2023-04-19)
--------------------------

- Fixed unicode error when reading fastq file

Version 1.0.1 (2023-03-28)
--------------------------

- Fixed invalid uppercase when iterating fastx

Version 1.0.0 (2023-03-24)
--------------------------

- Added support for fasta header without space
- Fixed some files missing in pypi tar.gz file

Version 0.9.1 (2022-12-31)
--------------------------

- Fixed unicode decode error when parsing large fasta/q file
- Fixed sequence retrival error when using sequence object from loop after break

Version 0.9.0 (2022-12-30)
--------------------------

- Added support for Python3.10, 3.11
- Added support for aarch64 and musllinux
- Added using tab as fasta sequence name splitter
- Fixed repeat sequence comment error
- Fixed the quality score parsing error from fastq
- Fixed the reference of sequence returned from function

Older versions
--------------

Version 0.8.4 (2021-06-30)
^^^^^^^^^^^^^^^^^^^^^^^^^^

- Added slice feature to FastaKeys
- Fixed FastaKeys and FastqKeys iteration memory leak
- Optimized FastaKeys and FastqKeys creation

Version 0.8.3 (2021-04-25)
^^^^^^^^^^^^^^^^^^^^^^^^^^

- Fixed Fastx iteration for next function
- Fixed Fastx uppercase for reading fasta

Version 0.8.2 (2021-01-02)
^^^^^^^^^^^^^^^^^^^^^^^^^^

- Fixed sample segfault error caused by fastq iteration error
- Fixed gzip index import error in multiple processes
- Fixed fastq iteration segfault error with full_name=True
- Fixed all objects iteration to support built-in next function

Version 0.8.1 (2020-12-16)
^^^^^^^^^^^^^^^^^^^^^^^^^^

- Fixed pip install error from source code
- Removed support for python39 32bit due to dll load error

Version 0.8.0 (2020-12-15)
^^^^^^^^^^^^^^^^^^^^^^^^^^

- Added Fastx object as a simple sequence iterator
- Added FastqKeys object to obtain read names
- Added full_name option to Fastq object
- Added support for Python 3.9
- Fixed Fasta object error identifier order
- Optimized speed of containing test and iteration
- Changed Identifier object to FastaKeys object

Version 0.7.0 (2020-09-20)
^^^^^^^^^^^^^^^^^^^^^^^^^^

- Added support for extracting flank sequences
- Added support for indexing super large gzip file
- Reduced memory consumption when building gzip index
- Improved the speed of random access to reads from fastq
- Fixed sequence dealloc error cuasing no fasta delloc trigger
- Fixed fastq max and min quality score return value

Version 0.6.17 (2020-08-31)
^^^^^^^^^^^^^^^^^^^^^^^^^^^

- Fixed gzip index loading error when no write permission

Version 0.6.16 (2020-08-27)
^^^^^^^^^^^^^^^^^^^^^^^^^^^

- Increased the buff size of kseq to speedup sequence iteration
- Removed warning message from fasta.c when building full index

Version 0.6.15 (2020-08-25)
^^^^^^^^^^^^^^^^^^^^^^^^^^^

- Fixed key_func error caused by free operation
- Fixed full name error when reading sequence without whitespace in names
- Fixed a hidden bug in fasta/q iteration when reading attributes (not seq)
- Fixed fasta/fastq size and sequence count error on Windows when parsing large file
- Fixed zlib 2gb and 4gb limit on windows x64 to support large file
- Reduced seek point span size to speedup random access from gzip file

Version 0.6.14 (2020-07-31)
^^^^^^^^^^^^^^^^^^^^^^^^^^^

- Added support for using full header as identifier without building index
- Improved the speed of fasta sequence iteration
- Improved the speed of gzipped fastq read iteration
- Fixed a bug in fastq read reader

Version 0.6.13 (2020-07-09)
^^^^^^^^^^^^^^^^^^^^^^^^^^^

- Fixed fastq read iteration error
- Fixed fastq cache buffer reader
- Added cache for mean, median and N50 length
- Speedup fasta iteration by reduced seeks

Version 0.6.12 (2020-06-14)
^^^^^^^^^^^^^^^^^^^^^^^^^^^

- Fixed DeprecationWarning on py38 caused by '#' formats args
- Fixed some memory leak bugs
- Cached sequence name to speedup fetch method
- Used random string as gzip index temp file to support multiple processes


Version 0.6.11 (2020-05-18)
^^^^^^^^^^^^^^^^^^^^^^^^^^^

- Fixed iteration error on Windows
- Fixed test error on Windows
- Fixed fastq composition error on 32bit OS
- Improved the speed of fasta identifier sort and filter

Version 0.6.10 (2020-04-22)
^^^^^^^^^^^^^^^^^^^^^^^^^^^

- Improved the speed of sequence reading 
- Improved the speed of sequence line iteration
- Added avglen, minlen, maxlen, minqual and maxqual to Fastq object
- Fixed read retrieval error
- Fixed some hidden memory leaks
- Changed fastq index file structure to save more information

Version 0.6.9 (2020-04-12)
^^^^^^^^^^^^^^^^^^^^^^^^^^

- Added buffreader to improve speed for reading from gzipped file
- Added extract subcommand to extract sequences from fasta/q file
- Added build subcommand to just build index
- Changed info subcommand output to a tab seperated table
- Changed Fastq object composition parameter to full_index

Version 0.6.8 (2020-03-14)
^^^^^^^^^^^^^^^^^^^^^^^^^^

- Fixed large offset seek error on windows
- Fixed PyUnicode_AsUTF8 const char type warning
- Changed sequence read line by line function
- Changed gzread to fread for fastq information

Version 0.6.7 (2020-03-03)
^^^^^^^^^^^^^^^^^^^^^^^^^^

- Added check for fasta/q format when open file
- Added benchmark scripts for evaluating performance
- Speed up the fasta/q object iteration
- Optimzed str length warning caused by strlen

Version 0.6.6 (2020-02-15)
^^^^^^^^^^^^^^^^^^^^^^^^^^

- Fixed incorrect sliced sequence name
- Fixed seq,identifier,read object memory dealloc
- Changed description text into description length in index file

Version 0.6.5 (2020-01-31)
^^^^^^^^^^^^^^^^^^^^^^^^^^

- Reduced memory usage when building index for large fasta
- Removed rebuild_index method from Fasta object due to segmentation fault
- Optimized compatibility between sqlite3 and python GIL

Version 0.6.4 (2020-01-14)
^^^^^^^^^^^^^^^^^^^^^^^^^^

- Fixed last sequence fetching error caused by missing \n
- Improved fasta/q object key error message to make it more human

Version 0.6.3 (2020-01-08)
^^^^^^^^^^^^^^^^^^^^^^^^^^

- Added .raw attribute to sequence object to get seq raw string
- Added .raw attribute to read object to get read raw string
- Added .description to read object to get full header line
- Added iteration for sequence object from FASTA object
- Added iteration for tuple from FASTQ object
- Changed FASTA class parameter composition to full_index

Version 0.6.2 (2020-01-04)
^^^^^^^^^^^^^^^^^^^^^^^^^^

- Fixed sample sequence index error
- Fixed ci deploy error

Version 0.6.1 (2020-01-03)
^^^^^^^^^^^^^^^^^^^^^^^^^^

- Added sample sequences command line
- Added get subsequence command line

Version 0.6.0 (2020-01-02)
^^^^^^^^^^^^^^^^^^^^^^^^^^

- Fixed FASTA object parameter error
- Fixed identifier sprintf warning
- Fixed fasta description end \r retained
- Fixed error byte length when slice sequence
- Removed support for python2.7 and python3.4
- Removed python2 compat
- Disabled export gzip index when building memory index

Version 0.5.10 (2019-11-20)
^^^^^^^^^^^^^^^^^^^^^^^^^^^

- Added identifier filter function
- Remove tp_new for Read, Sequence and Identifier
- Fixed module method error

Version 0.5.9 (2019-11-17)
^^^^^^^^^^^^^^^^^^^^^^^^^^

- Added get longest and shortest sequence object
- Added composition argument to speedup getting GC content
- Added memory index to keep index in memory rather than local file
- Fixed command line error
- Changed sqlite to higher version
- Removed journal_mode OFF
- Speedup index building

Version 0.5.8 (2019-11-10)
^^^^^^^^^^^^^^^^^^^^^^^^^^

- Fixed fasta NL function parameter check
- Fixed read id error when fastq iteration

Version 0.5.7 (2019-11-09)
^^^^^^^^^^^^^^^^^^^^^^^^^^

- Fixed SystemError caused caused by Python 2.7 seperated int and long type
- Fixed String type check on Python 2.7
- Fixed objects memory deallocation

Version 0.5.6 (2019-11-08)
^^^^^^^^^^^^^^^^^^^^^^^^^^

- Optimized random access from plain file
- Reduced memory consumption

Version 0.5.5 (2019-11-07)
^^^^^^^^^^^^^^^^^^^^^^^^^^

- Added Support for IUPAC code complement
- Speedup reverse complement
- Speedup space removing and uppercase


Version 0.5.4 (2019-11-04)
^^^^^^^^^^^^^^^^^^^^^^^^^^

- Added guess fasta type (DNA, RNA, protein)
- Added support for calculating protein sequence composition
- Optimized the speed of index building
- Calculate sequence composition when get gc content or composition
- Fixed char return in python 2.7

Version 0.5.3 (2019-10-23)
^^^^^^^^^^^^^^^^^^^^^^^^^^

- Added support for coverting fastq to fasta
- Updated command line interface docs
- Fixed command line entry points

Version 0.5.2 (2019-10-18)
^^^^^^^^^^^^^^^^^^^^^^^^^^

- Fixed command line interface running error

Version 0.5.1 (2019-10-17)
^^^^^^^^^^^^^^^^^^^^^^^^^^

- Added key function for custom sequence identifier
- Optimized speed of fasta indexing
- Fixed bool args parsing error in py2.7

Version 0.5.0 (2019-10-13)
^^^^^^^^^^^^^^^^^^^^^^^^^^

- Added support for python 2.7 and 3.4
- Added command line tool to manipulate fasta and fastq file
- Added gzip attribute to fasta and fastq object to check whether compressed
- Added sort function for identifier object
- Fixed python bool argument parsing error caused by uint16_t
- Fixed identifier sort key initialization

Version 0.4.1 (2019-10-05)
^^^^^^^^^^^^^^^^^^^^^^^^^^

- Fixed fastq quality encoding system guesser
- Fixed gzip index insertion error

Version 0.4.0 (2019-09-29)
^^^^^^^^^^^^^^^^^^^^^^^^^^

- Added support for parsing FASTQ
- Added random access to reads from FASTQ

Version 0.3.10 (2019-09-27)
^^^^^^^^^^^^^^^^^^^^^^^^^^^

- Fixed GC skew exception caused by mixing unsigned with signed for division  

Version 0.3.9 (2019-09-26)
^^^^^^^^^^^^^^^^^^^^^^^^^^

- Fixed sequence read line by line error
- Fixed last sequence build index error when fasta file ended without \n
- Fixed GC skew error

Version 0.3.8 (2019-09-25)
^^^^^^^^^^^^^^^^^^^^^^^^^^

- Fixed large offset became negative error
- Fixed slice step
- Fixed uncorrect median length
- Fixed strand compare error
- Added GC skew calculation
- Updated test script

Version 0.3.7 (2019-09-24)
^^^^^^^^^^^^^^^^^^^^^^^^^^

- Changed int type to standard type
- Added support for processing large fasta file
- Added id number for each sequence
- Fixed SQL fetch error
- Used 50 as default value of nl to calculate N50 and L50

Version 0.3.6 (2019-09-20)
^^^^^^^^^^^^^^^^^^^^^^^^^^

- Added support for searching subsequence from a sequence
- Added support for checking subsequence weather in a sequence
- Fixed gzip index import error
- Fixed subsequence parent length for full sequence extraction

Version 0.3.5 (2019-09-08)
^^^^^^^^^^^^^^^^^^^^^^^^^^

- Fixed unicode error caused by sqlite3_finalize 

Version 0.3.4 (2019-09-07)
^^^^^^^^^^^^^^^^^^^^^^^^^^

- Fixed seq description unicode string error

Version 0.3.3 (2019-09-07)
^^^^^^^^^^^^^^^^^^^^^^^^^^

- Fixed sequence description encoding error
 
Version 0.3.2 (2019-09-07)
^^^^^^^^^^^^^^^^^^^^^^^^^^

Deleted

Version 0.3.1 (2019-09-07)
^^^^^^^^^^^^^^^^^^^^^^^^^^

- Added support for geting sequence description

Version 0.3.0 (2019-09-07)
^^^^^^^^^^^^^^^^^^^^^^^^^^

- Added read sequence from fasta file line by line
- Added support for calculating assembly N50 and L50
- Added support for calculating median and average length
- Added support for getting longest and shortest sequence
- Added support for calculating counts of sequence
- removed support for Python34

Version 0.2.11 (2019-08-31)
^^^^^^^^^^^^^^^^^^^^^^^^^^^

- Support for Python 3.4

Version 0.2.10 (2019-08-28)
^^^^^^^^^^^^^^^^^^^^^^^^^^^

- Changed fseek and fread into gzseek and gzread
- Fixed sequence cache name comparision
- Fixed last sequence read error without line end
- Fixed subsequence slice error in normal FASTA file

Version 0.2.9 (2019-08-27)
^^^^^^^^^^^^^^^^^^^^^^^^^^

- Fixed bad line calculate error
- Changed rewind to fseek for subsequence extraction

Version 0.2.8 (2019-08-26)
^^^^^^^^^^^^^^^^^^^^^^^^^^

- Changed kseq.h library from li to attractivechaos
- Improved fasta parser

Version 0.2.7 (2019-08-26)
^^^^^^^^^^^^^^^^^^^^^^^^^^

- Fixed no gzip index wrote to sqlite index file

Version 0.2.6 (2019-08-26)
^^^^^^^^^^^^^^^^^^^^^^^^^^

- Optimized speed of gzip random access

Version 0.2.5 (2019-08-25)
^^^^^^^^^^^^^^^^^^^^^^^^^^

- Fixed segmentation fault raised when loading gzip index
- Changed fasta object method get_seq to fetch

Version 0.2.4 (2019-08-25)
^^^^^^^^^^^^^^^^^^^^^^^^^^

- Fixed fasta iter error after building new index

Version 0.2.3 (2019-08-24)
^^^^^^^^^^^^^^^^^^^^^^^^^^

- Fixed fasta iter error when end of file is not \n

Version 0.2.2 (2019-07-19)
^^^^^^^^^^^^^^^^^^^^^^^^^^

- Fixed identifier contain error

Version 0.2.1 (2019-07-15)
^^^^^^^^^^^^^^^^^^^^^^^^^^

- Fixed sequence name always end with 0
- Fixed fasta iterable for flat fasta

Version 0.2.0 (2019-07-09)
^^^^^^^^^^^^^^^^^^^^^^^^^^

- First release to PyPI
