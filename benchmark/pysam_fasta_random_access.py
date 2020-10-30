import sys
import pysam

idfile, fafile = sys.argv[1:]

fa = pysam.FastaFile(fafile)

with open(idfile) as fh:
	for line in fh:
		seqid = line.strip()
		s = str(fa[seqid])
		print(s)
