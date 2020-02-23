import sys
import pyfaidx

print(sys.argv)

chrom, fafile = sys.argv[1:]

fa = pyfaidx.Fasta(fafile)

subseq = str(-fa[chrom][:])

print(">{}/rc pyfastx\n{}\n".format(chrom, subseq))
