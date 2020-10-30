import sys
from Bio import SeqIO
fq = SeqIO.index_db("{}.db".format(sys.argv[2]))

with open(sys.argv[1]) as fh:
	for line in fh:
		name = line.strip()
		print(fq[name].seq)

