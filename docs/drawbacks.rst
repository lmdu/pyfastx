Drawbacks
=========

If you intensively check sequence names exists in FASTA file using ``in`` operator on FASTA object like:

.. code:: python

	>>> fa = pyfastx.Fasta('tests/data/test.fa.gz')
	>>> # Suppose seqnames has 100000 names
	>>> for seqname in seqnames:
	>>>     if seqname in fa:
	>>>	        do something

This will take a long time to finish. Because, pyfastx does not load the index into memory, the ``in`` operating is corresponding to sql query existence from index database. The faster alternative way to do this is:

.. code:: python

	>>> fa = pyfastx.Fasta('tests/data/test.fa.gz')
	>>> # load all sequence names into a set object
	>>> all_names = set(fa.keys())
	>>> for seqname in seqnames:
	>>>     if seqname in all_names:
	>>>	        do something
