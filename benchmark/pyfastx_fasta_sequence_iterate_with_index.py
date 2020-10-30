import sys
import pyfastx

for seq in pyfastx.Fasta(sys.argv[1]):
	_ = seq.seq
