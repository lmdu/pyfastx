FASTA
=====

Read FASTA file
---------------

The fastest way to parse flat or gzipped FASTA file without building index.

.. code:: python

	>>> import pyfastx
	>>> for name, seq in pyfastx.Fasta('test/data/test.fa.gz', build_index=False):
	>>> 	print(name, seq)

Read flat or gzipped FASTA file and build index, support for random access to FASTA.

.. code:: python

	>>> import pyfastx
	>>> fa = pyfastx.Fasta('test/data/test.fa.gz')
	>>> fa
	<Fasta> test/data/test.fa.gz contains 211 seqs

.. note::
	Building index may take some times. The time required to build index depends on the size of FASTA file. If index built, you can randomly access to any sequences in FASTA file.

Get FASTA information
---------------------

.. code:: python

	>>> # get sequence counts in FASTA
	>>> len(fa)
	211

	>>> # get total sequnce length of FASTA
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
	{'A': 24534, 'C': 18694, 'G': 18855, 'T': 24179, 'N': 0}

Get longest and shortest sequence
---------------------------------

New in ``pyfastx`` 0.3.0

.. code:: python

	>>> # get longest sequence (name, length)
	>>> fa.longest
	('JZ822609.1', 821)

	>>> # get shortest sequence (name, length)
	>>> fa.shortest
	('JZ822617.1', 118)

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

The fastest way to parse plain or gzipped FASTQ file without building index.

.. code:: python

	>>> import pyfastx
	>>> for read in pyfastx.Fastq('tests/data/test.fq.gz', build_index=False):
	>>> 	print(read.name, read.seq, read.qual)

Read plain or gzipped file and build index, support for random access to reads from FASTQ.

.. code:: python

	>>> import pyfastx
	>>> fq = pyfastx.Fastq('tests/data/test.fq.gz')
	>>> fq
	<Fastq> tests/data/test.fq.gz contains 100 reads

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
	>>> fq.guess
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

	>>> # get read length
	>>> len(r)
	150

	>>> # get read sequence
	>>> r.seq
	'CGAGGAAATCGACGTCACCGATCTGGAAGCCCTGCGCGCCCATCTCAACCAGAAATGGGGTGGCCAGCGCGGCAAGCTGACCCTGCTGCCGTTCCTGGTCCGCGCCATGGTCGTGGCGCTGCGCGACTTCCCGCAGTTGAACGCGCGCTA'

	>>> # get read quality ascii string
	>>> r.qual
	'FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FF,FFFFFFFFFFFFFFFFFFFFFFFFFF,F:FFFFFFFFF:'

	>>> # get read quality integer value, ascii - 33 or 64
	>>> r.quali
	[37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 25, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 25, 37, 37, 11, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 11, 37, 25, 37, 37, 37, 37, 37, 37, 37, 37, 37, 25]

	>>> # get read length
	>>> len(r)
	150

Identifier
==========

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
	>>> 	print(name)

	>>> # convert to a list
	>>> list(ids)

Sort identifiers
----------------

Sort identifiers by id, name, or length for iteration

New in ``pyfastx`` 0.5.0

.. code:: python

	>>> # sort identifiers by length with descending order 
	>>> for name in ids.sort(key='length', reverse=True):
	>>> 	print(name)

	>>> # sort identifiers by name with ascending order
	>>> for name in ids.sort(key='name'):
	>>> 	print(name)

	>>> # sort identifiers by id with descending order
	>>> for name in ids.sort(key='id', reverse=True)
	>>> 	print(name)
