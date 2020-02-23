import sys
import pyfastx

chrom, fafile = sys.argv[1:]

fa = pyfastx.Fasta(fafile)

subseq = str(fa[chrom][99:200])

print(">{}:100-200 pyfastx\n{}\n".format(chrom, subseq))
