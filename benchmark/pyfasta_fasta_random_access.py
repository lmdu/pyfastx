import sys
import pyfasta

idfile, fafile = sys.argv[1:]

fa = pyfasta.Fasta(fafile)

with open(idfile) as fh:
	for line in fh:
		seqid = line.strip()
		s = str(fa[seqid])
		print(s)
