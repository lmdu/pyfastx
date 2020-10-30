import sys
import pyfaidx

idfile, fafile = sys.argv[1:]

fa = pyfaidx.Fasta(fafile)

with open(idfile) as fh:
	for line in fh:
		seqid, start, end = line.strip().split()
		s = str(fa[seqid][int(start):int(end)])
		print(s)
