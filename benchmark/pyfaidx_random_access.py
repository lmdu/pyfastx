import sys
import pyfaidx

chrom, fafile = sys.argv[1:]

fa = pyfaidx.Fasta(fafile)

subseq = str(fa[chrom][99:200])

print(">{}:100-200 pyfaidx\n{}\n".format(chrom, subseq))
