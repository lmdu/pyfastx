import sys
import pyfasta

chrom, fafile = sys.argv[1:]

fa = pyfasta.Fasta(fafile, key_fn=lambda key: key.split()[0])

subseq = str(fa[chrom][99:200])

print(">{}:100-200 pyfasta\n{}\n".format(chrom, subseq))
