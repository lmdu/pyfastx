API Reference
=============

pyfastx.version
---------------

.. py:function:: pyfastx.version(debug=False)

	Get current version of pyfastx

	:param bool debug: if true, return versions of pyfastx, zlib, sqlite3 and zran.

	:return: version of pyfastx

	:rtype: str

.. py:function:: pyfastx.gzip_check(file_name)

	New in pyfastx 0.5.4

	Check file is gzip compressed or not

	:param str file_name: the path of input file

	:return: Ture if file is gzip compressed else False

	:rtype: bool

.. py:function:: pyfastx.reverse_complement(seq)

	New in pyfastx 2.0.0

	get reverse complement sequence of given DNA sequence

	:param str seq: DNA sequence

	:return: reverse complement sequence

	:rtype: str

pyfastx.Fasta
-------------

.. py:class:: pyfastx.Fasta(file_name, index_file=None, uppercase=True, build_index=True, full_index=False, full_name=False, memory_index=False, key_func=None)

	Read and parse fasta files. Fasta can be used as dict or list, you can use index or sequence name to get a sequence object, e.g. ``fasta[0]``, ``fasta['seq1']``

	:param str file_name: the file path of input FASTA file

	:param str index_file: the index file of FASTA file, default using index file with extension of .fxi in the same directory of FASTA file, New in 2.0.0

	:param bool uppercase: always output uppercase sequence, default: ``True``

	:param bool build_index: build index for random access to FASTA sequence, default: ``True``. If build_index is False, iteration will return a tuple (name, seq); If build_index is True, iteration will return a sequence object.

	:param bool full_index: calculate character (e.g. A, T, G, C) composition when building index, this will improve the speed of GC content extracting. However, it will take more time to build index, default: ``False``

	:param bool full_name: use the full header line instead of the part before first whitespace as the identifier of sequence, even in mode without building index. New in 0.6.14, default: ``False``

	:param bool memory_index: if memory_index is True, the fasta index will be kept in memory and do not generate a index file, default: ``False``

	:param function key_func: new in 0.5.1, key function is generally a lambda expression to split header and obtain a shortened identifer, default: ``None``

	:return: Fasta object

	.. py:attribute:: file_name

		FASTA file path

	.. py:attribute:: size

		total length of sequences in FASTA file

	.. py:attribute:: type

		New in ``pyfastx`` 0.5.4

		get fasta type, return DNA, RNA, protein, or unknown

	.. py:attribute:: is_gzip

		New in pyfastx 0.5.0

		return True if fasta is gzip compressed else return False

	.. py:attribute:: gc_content

		GC content of whole sequences in FASTA file, return a float value

	.. py:attribute:: gc_skew

		GC skew of whole sequences in FASTA file, learn more about `GC skew <https://en.wikipedia.org/wiki/GC_skew>`_

		New in ``pyfastx`` 0.3.8

	.. py:attribute:: composition

		nucleotide composition in FASTA file, a dict contains counts of A, T, G, C and N (unkown nucleotide base)

	.. py:attribute:: longest

		get longest sequence in FASTA file, return a Sequence object

		New in ``pyfastx`` 0.3.0

	.. py:attribute:: shortest

		get shortest sequence in FASTA file, return a Sequence object

		New in ``pyfastx`` 0.3.0

	.. py:attribute:: mean

		get average length of sequences in FASTA file

		New in ``pyfastx`` 0.3.0

	.. py:attribute:: median

		get median length of sequences in FASTA file

		New in ``pyfastx`` 0.3.0

	.. py:method:: fetch(chrom, intervals, strand='+')

		truncate subsequences from a given sequence by a start and end coordinate or a list of coordinates. This function will cache the full sequence into memory, and is suitable for extracting large numbers of subsequences from specified sequence.

		:param str chrom: chromosome name or sequence name

		:param list/tuple intervals: list of [start, end] coordinates

		:param str strand: sequence strand, ``+`` indicates sense strand, ``-`` indicates antisense strand, default: '+'

		.. note::

			intervals can be a list or tuple with start and end position e.g. (10, 20).
			intervals also can be a list or tuple with multiple coordinates e.g. [(10, 20), (50,70)]

		:return: sliced subsequences

		:rtype: str

	.. py:method:: flank(chrom, start, end, flank_length=50, use_cache=False)

		Get the flank sequence of given subsequence with start and end. New in 0.7.0

		:param str chrom: chromosome name or sequence name

		:param int start: 1-based start position of subsequence on chrom

		:param int end: 1-based end position of subsequence on chrom

		:param int flank_length: length of flank sequence, default 50

		:param bool use_cache: cache the whole sequence

		.. note::

			If you want to extract flank sequence for large numbers of subsequences from the same sequence. Use ``use_cache=True`` will greatly improve the speed

		:return: left flank and right flank sequence

		:rtype: tuple

	.. py:method:: build_index()

		build index for FASTA file

	.. py:method:: keys()

		get all names of sequences

		:return: an FastaKeys object

	.. py:method:: count(n)

		get counts of sequences whose length >= n bp

		New in ``pyfastx`` 0.3.0

		:param int n: number of bases

		:return: sequence counts

		:rtype: int

	.. py:method:: nl(quantile)

		calculate assembly N50 and L50, learn more about `N50,L50 <https://www.molecularecologist.com/2017/03/whats-n50/>`_

		New in ``pyfastx`` 0.3.0

		:param int quantile: a number between 0 and 100, default 50

		:return: (N50, L50)

		:rtype: tuple

pyfastx.Sequence
----------------

.. py:class:: pyfastx.Sequence

	Readonly sequence object generated by fasta object, Sequence can be treated as a list and support slicing e.g. ``seq[10:20]``

	.. py:attribute:: id

		sequence id or order number in FASTA file

	.. py:attribute:: name

		sequence name

	.. py:attribute:: description

		Get sequence description after name in sequence header

		New in ``pyfastx`` 0.3.1

	.. py:attribute:: start

		start position of sequence

	.. py:attribute:: end

		end position of sequence

	.. py:attribute:: gc_content

		GC content of current sequence, return a float value

	.. py:attribute:: gc_skew

		GC skew of current sequence, learn more about `GC skew <https://en.wikipedia.org/wiki/GC_skew>`_

	.. py:attribute:: composition

		nucleotide composition of sequence, a dict contains counts of A, T, G, C and N (unkown nucleotide base)

	.. py:attribute:: raw

		get the raw string (with header line and sequence lines) of sequence as it appeared in file

		New in ``pyfastx`` 0.6.3

	.. py:attribute:: seq

		get the string of sequence in sense strand

	.. py:attribute:: reverse

		get the string of reversed sequence

	.. py:attribute:: complement

		get the string of complement sequence

	.. py:attribute:: antisense

		get the string of sequence in antisense strand, corresponding to reversed and complement sequence

	.. py:method:: search(subseq, strand='+')

		Search for subsequence from given sequence and get the start position of the first occurrence

		New in ``pyfastx`` 0.3.6

		:param str subseq: a subsequence for search

		:param str strand: sequence strand + or -, default +

		:return: if found subsequence return one-based start position, if not return None

		:rtype: int or None

pyfastx.Fastq
-------------

New in ``pyfastx`` 0.4.0

.. py:class:: pyfastx.Fastq(file_name, index_file=None, phred=0, build_index=True, full_index=False)

	Read and parse fastq file

	:param str file_name: input FASTQ file path

	:param str index_file: the index file of FASTQ file, default using the index file with extension of .fxi in the same directory of FASTQ file. New in 2.0.0

	:param bool build_index: build index for random access to FASTQ reads, default: ``True``. If build_index is False, iteration will return a tuple (name, seq, qual); If build_index is True, iteration will return a read object

	:param bool full_index: calculate character (e.g. A, T, G, C) composition when building index, this will improve the speed of GC content extracting. However, it will take more time to build index, default: ``False``

	:param int phred: phred was used to convert quality ascii to quality int value, usually is 33 or 64, default ``33``

	:return: Fastq object

	.. py:attribute:: file_name

		FASTQ file path

	.. py:attribute:: size

		total bases in FASTQ file

	.. py:attribute:: is_gzip

		New in pyfastx 0.5.0

		return True if fasta is gzip compressed else return False

	.. py:attribute:: gc_content

		GC content of whole FASTQ file

	.. py:attribute:: avglen

		New in ``pyfastx`` 0.6.10

		get average length of reads

	.. py:attribute:: maxlen

		New in ``pyfastx`` 0.6.10

		get maximum length of reads

	.. py:attribute:: minlen

		New in ``pyfastx`` 0.6.10

		get minimum length of reads

	.. py:attribute:: maxqual

		New in ``pyfastx`` 0.6.10

		get maximum quality score of bases

	.. py:attribute:: minqual

		New in ``pyfastx`` 0.6.10

		get minimum quality score of bases

	.. py:attribute:: composition

		base composition in FASTQ file, a dict contains counts of A, T, G, C and N (unkown nucleotide base)

	.. py:attribute:: phred

		get phred value

	.. py:attribute:: encoding_type

		New in ``pyfastx`` 0.4.1

		Guess the quality encoding type used by FASTQ sequence file

	.. py:method:: build_index()

		Build index for fastq file when build_index set to False

	.. py:method:: keys()

		New in ``pyfastx`` 0.8.0

		Get all the names of reads in fastq file

		:return: an FastqKeys object

pyfastx.Read
------------

New in ``pyfastx`` 0.4.0

.. py:class:: pyfastx.Read

	Readonly read object for obtaining read information, generated by fastq object

	.. py:attribute:: id

		read id or order number in FASTQ file

	.. py:attribute:: name

		read name excluding '@'

	.. py:attribute:: description

		get the full header line of read

	.. py:attribute:: raw

		get the raw string (with header, sequence, comment and quality lines) of read as it appeared in file

		New in ``pyfastx`` 0.6.3

	.. py:attribute:: seq

		get read sequence string

	.. py:attribute:: qual

		get read quality ascii string

	.. py:attribute:: quali

		get read quality integer value (ascii - phred), return a list

pyfastx.Fastx
-------------

.. py:class:: pyfastx.Fastx(file_name, format="auto", uppercase=False)

	New in ``pyfastx`` 0.8.0. A python binding of kseq.h, provide a simple api for iterating over sequences in fasta/q file

	:param str file_name: input fasta or fastq file path

	:param str format: the input file format, can be "fasta" or "fastq", default: "auto", automatically detect the format of sequence file

	:param bool uppercase: always output uppercase sequence, only work for fasta file, default: False

	:return: Fastx object

pyfastx.FastaKeys
------------------

.. py:class:: pyfastx.FastaKeys

	FastaKeys is a readonly and list-like object, contains all names of sequences

	.. py:method:: sort(by="id", reverse=False)

		Sort keys by sequence id, name or length for iteration

		New in ``pyfastx`` 0.5.0

		:param str by: order by id, name, or length, default is id

		:param bool reverse: used to flag descending sorts, default is False

		:return: FastaKeys object itself

	.. py:method:: filter(*filters)

		Filter keys by sequence name and length for iteration

		:param list filters: filters generated by comparison like ids > 500 or ids % 'seq1', where ids is a Identifier object

		:return: FastaKeys object itself

	.. py:method:: reset()

		Clear all filters and sort order

		:return: FastaKeys object itself

pyfastx.FastqKeys
------------------

.. py:class:: pyfastx.FastqKeys

	New in ``pyfastx`` 0.8.0. FastqKeys is a readonly and list-like object, contains all names of reads
