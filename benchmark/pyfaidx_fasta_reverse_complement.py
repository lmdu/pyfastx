import sys
import pyfaidx

chrom, fafile = sys.argv[1:]

fa = pyfaidx.Fasta(fafile)
seq = -fa[chrom][:]
print(seq.seq)
