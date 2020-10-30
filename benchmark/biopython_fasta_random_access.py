import sys
from Bio import SeqIO

idfile, fafile = sys.argv[1:]

fa = SeqIO.index_db("{}.db".format(fafile))

with open(idfile) as fh:
	for line in fh:
		seqid = line.strip()
		s = fa[seqid].seq
		print(s)
