import sys
import pysam

chrom, fafile = sys.argv[1:]

fa = pysam.FastaFile(fafile)

subseq = fa.fetch(chrom, 100, 200)

print(">{}:100-200 pysam\n{}\n".format(chrom, subseq))