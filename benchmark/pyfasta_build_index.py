import sys
import pyfasta

pyfasta.Fasta(sys.argv[1], key_fn=lambda key: key.split()[0])
