import sys
import pyfaidx

idfile, fafile = sys.argv[1:]

fa = pyfaidx.Fasta(fafile)

with open(idfile) as fh:
	for line in fh:
		seqid = line.strip()
		s = str(fa[seqid])
		print(s)
