import sys
import pyfastx

for name, seq in pyfastx.Fasta(sys.argv[1], build_index=False):
	pass
