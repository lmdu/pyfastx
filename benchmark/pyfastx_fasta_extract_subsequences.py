import sys
import pyfastx

idfile, fafile = sys.argv[1:]

fa = pyfastx.Fasta(fafile)

with open(idfile) as fh:
	for line in fh:
		seqid, start, end = line.strip().split()
		s = fa[seqid][int(start):int(end)].seq
		print(s)
