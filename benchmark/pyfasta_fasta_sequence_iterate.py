import sys
import pyfasta

fa = pyfasta.Fasta(sys.argv[1])
for name in fa:
	_ = str(fa[name])
