import sys
from Bio import SeqIO

idfile, fafile = sys.argv[1:]

fa = SeqIO.index_db("{}.db".format(fafile))

with open(idfile) as fh:
	for line in fh:
		seqid, start, end = line.strip().split()
		s = fa[seqid][int(start):int(end)].seq
		print(s)
