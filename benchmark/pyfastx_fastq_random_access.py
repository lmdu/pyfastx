import sys
import pyfastx
fq = pyfastx.Fastq(sys.argv[2])

with open(sys.argv[1]) as fh:
	for line in fh:
		name = line.strip()
		print(fq[name].seq)
