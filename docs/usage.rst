FASTX
=====

New in ``pyfastx`` 0.8.0.

Iterate over sequences in FASTA
-------------------------------

When iterating over sequences on FASTX object, a tuple ``(name, seq)`` will be returned.

.. code:: python

	>>> fa = pyfastx.Fastx('tests/data/test.fa')
	>>> for name,seq in fa:
	>>> 	print(name)
	>>> 	print(seq)

	>>> #always output uppercase sequence
	>>> for item in pyfastx.Fastx('tests/data/test.fa', uppercase=True):
	>>> 	print(item)

	>>> #Manually specify sequence format
	>>> for item in pyfastx.Fastx('tests/data/test.fa', format="fasta"):
	>>> 	print(item)

If you want the sequence comment, you can set comment to True, New in ``pyfastx`` 0.9.0.

.. code:: python

    >>> fa = pyfastx.Fastx('tests/data/test.fa.gz', comment=True)
    >>> for name,seq,comment in fa:
    >>>     print(name)
    >>>     print(seq)
    >>>     print(comment)

.. note::

	The comment is the content of header line after the first white space or tab character.

Iterate over reads in FASTQ
---------------------------

When iterating over reads on FASTX object, a tuple ``(name, seq, qual)`` will be returned.

.. code:: python

	>>> fq = pyfastx.Fastx('tests/data/test.fq')
	>>> for name,seq,qual in fq:
	>>> 	print(name)
	>>> 	print(seq)
	>>> 	print(qual)

If you want the read comment, you can set comment to True, New in ``pyfastx`` 0.9.0.

.. code:: python

    >>> fq = pyfastx.Fastx('tests/data/test.fq.gz', comment=True)
    >>> for name,seq,qual,comment in fq:
    >>>     print(name)
    >>>     print(seq)
    >>>     print(qual)
    >>>     print(comment)

.. note::

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
	Building index may take some time. The time required to build index depends on the size of FASTA file. If index built, you can randomly access to any sequences in FASTA file. The index file can be reused to save time when you read sequences from FASTA file next time.

FASTA records iteration
-----------------------

The fastest way to iterate plain or gzipped FASTA file without building index, the iteration will return a tuple contains name and sequence.

.. code:: python

	>>> import pyfastx
	>>> for name, seq in pyfastx.Fasta('test/data/test.fa.gz', build_index=False):
	>>> 	print(name, seq)

If you want to use full header line as sequence identifier without building index, you can do like this:

.. code:: python

	>>> import pyfastx
	>>> for name, seq in pyfastx.Fasta('test/data/test.fa', build_index=False, full_name=True):
	>>> 	print(name, seq)

You can also iterate sequence object from FASTA object like this:

.. code:: python

	>>> import pyfastx
	>>> for seq in pyfastx.Fasta('test/data/test.fa.gz'):
	>>> 	print(seq.name)
	>>> 	print(seq.seq)
	>>> 	print(seq.description)

Iteration with ``build_index=True`` (default) return sequence object which allows you to access attributes of sequence. New in pyfastx 0.6.3.

Get FASTA information
---------------------

.. code:: python

	>>> # get sequence counts in FASTA
	>>> len(fa)
	211

	>>> # get total sequence length of FASTA
	>>> fa.size
	86262

	>>> # get GC content of DNA sequences in FASTA
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
	>>> # New in pyfastx 0.5.4
	>>> fa.type
	'DNA'

	>>> # check fasta file is gzip compressed
	>>> # New in pyfastx 0.5.4
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

	>>> # get sequence median length
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

Get flank sequences
-------------------

New in ``pyfastx`` 0.7.0

Get flank sequences for the given subsequence, return a tuple with left and right flank sequence

.. code:: python

	>>> # get flank sequences with length of 20 for subsequence JZ822577:50-60
	>>> fa.flank('JZ822577.1', 50, 60, 20)
	('TCACTCAGGCTCTTTGTCAT', 'TAGGATATCGAGTATTCAAG')

	>>> # get flank sequences for a single base or SNP at position 100
	>>> fa.flank('JZ822577.1', 100, 100, 20)
	('GCTCATCGCTTTTGGTAATC', 'TTGCGGTGCATGCCTTTGCA')

	>>> # get flank sequences by buffer cache
	>>> fa.flank('JZ822577.1', 70, 90, flank_length=20, use_cache=True)
	('TTTAGTTTGACTAGGATATC', 'TTGGTAATCTTTGCGGTGCA')

.. note::

	The start and end position of subsequence were 1-based. When extracting flank for large numbers of subsequences from the same sequence, ``use_cache=True`` was recommended to improve speed.

Key function
------------

New in ``pyfastx`` 0.5.1

Sometimes your fasta will have a long header which contains multiple identifiers and description, for example:

``>JZ822577.1 contig1 cDNA library of flower petals in tree peony by suppression subtractive hybridization Paeonia suffruticosa cDNA, mRNA sequence``

In this case, either "JZ822577.1" or "contig1" could be used as the identifier.
You can specify the key function to select one as identifier.

.. code:: python

	>>> #default use JZ822577.1 as identifier
	>>> #specify key_func to select contig1 as identifer
	>>> fa = pyfastx.Fasta('tests/data/test.fa.gz', key_func=lambda x: x.split()[1])
	>>> fa
	<Fasta> tests/data/test.fa.gz contains 211 seqs

.. note::
	If the index file already existed, you should delete the previous index file, and then use key_func to create a new index file

Sequence
========

Get a sequence from FASTA
-------------------------

.. code:: python

	>>> # get sequence like dictionary
	>>> s1 = fa['JZ822577.1']
	>>> s1
	<Sequence> JZ822577.1 with length of 333

	>>> # get sequence like list
	>>> s2 = fa[2]
	>>> s2
	<Sequence> JZ822579.1 with length of 176

	>>> # get last sequence
	>>> s3 = fa[-1]
	>>> s3
	<Sequence> JZ840318.1 with length of 134

	>>> # check name weather in FASTA file
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

	>>> # get sequence description, New in pyfastx 0.3.1
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
	>>> 	print(name)
	>>> 	print(seq)
	>>> 	print(qual)

If you want to use full header line as read identifier without building index, you can do like this:

New in ``pyfastx`` 0.8.0

.. code:: python

	>>> import pyfastx
	>>> for name,seq,qual in pyfastx.Fastq('test/data/test.fq', build_index=False, full_name=True):
	>>> 	print(name, seq, qual)

You can also iterate read object from FASTQ object like this:

.. code:: python

	>>> import pyfastx
	>>> for read in pyfastx.Fastq('test/data/test.fq.gz'):
	>>> 	print(read.name)
	>>> 	print(read.seq)
	>>> 	print(read.qual)
	>>> 	print(read.quali)

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

	>>> # get the average length of reads
	>>> fq.avglen
	150.0

	>>> # get the maximum length of reads
	>>> fq.maxlen
	150

	>>> # get the minimum length of reads
	>>> fq.minlen
	150

	>>> # get the maximum quality score of bases
	>>> fq.maxqual
	70

	>>> # get the minimum quality score of bases
	>>> fq.minqual
	35

	>>> # Guess fastq quality encoding system
	>>> # New in pyfastx 0.4.1
	>>> fq.encoding_type
	['Sanger Phred+33', 'Illumina 1.8+ Phred+33']

Read
=====

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

FastaKeys
=========

Get fasta keys
---------------

Get all names of sequence as a list-like object

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
	>>> 	print(name)

	>>> # convert to a list
	>>> list(ids)

Sort keys
----------------

Sort keys by sequence id, name, or length for iteration

New in ``pyfastx`` 0.5.0

.. code:: python

	>>> # sort keys by length with descending order
	>>> for name in ids.sort(by='length', reverse=True):
	>>> 	print(name)

	>>> # sort keys by name with ascending order
	>>> for name in ids.sort(by='name'):
	>>> 	print(name)

	>>> # sort keys by id with descending order
	>>> for name in ids.sort(by='id', reverse=True)
	>>> 	print(name)

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

Clear filter and sort order
---------------------------

.. code:: python

	>>> # clear sort order and filters
	>>> ids.reset()
	<Identifier> contains 211 identifiers

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
