import sys
from Bio.SeqIO.QualityIO import FastqGeneralIterator

with open(sys.argv[1], 'rU') as handle:
	for name, seq, qual in FastqGeneralIterator(handle):
		print(seq)
