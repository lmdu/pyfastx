import sys
import pyfastx

fa = pyfastx.Fasta(sys.argv[1])
long_seq = fa.longest
print(fa.size, len(fa), long_seq.name, len(long_seq))
