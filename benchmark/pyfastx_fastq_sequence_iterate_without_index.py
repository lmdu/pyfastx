import sys
import pyfastx

for name,seq,qual in pyfastx.Fastq(sys.argv[1], build_index=False):
	print(seq)
