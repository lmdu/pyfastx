import sys
import pyfastx

idfile, fafile = sys.argv[1:]

fa = pyfastx.Fasta(fafile)

with open(idfile) as fh:
	for line in fh:
		seqid = line.strip()
		s = fa[seqid].seq
		print(s)
