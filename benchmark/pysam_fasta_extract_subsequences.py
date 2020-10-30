import sys
import pysam

idfile, fafile = sys.argv[1:]

fa = pysam.FastaFile(fafile)

with open(idfile) as fh:
	for line in fh:
		seqid, start, end = line.strip().split()
		s = fa.fetch(seqid, int(start), int(end))
		print(s)
