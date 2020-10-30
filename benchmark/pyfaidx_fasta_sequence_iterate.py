import sys
import pyfaidx

for seq in pyfaidx.Fasta(sys.argv[1]):
	print(str(seq))
