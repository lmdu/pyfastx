import sys
import pyfasta

idfile, fafile = sys.argv[1:]

fa = pyfasta.Fasta(fafile)

with open(idfile) as fh:
	for line in fh:
		seqid, start, end = line.strip().split()
		s = str(fa[seqid][int(start):int(end)])
		print(s)
