import sys
import pyfastx

chrom, fafile = sys.argv[1:]

fa = pyfastx.Fasta(fafile)
subseq = fa[chrom].seq
print(len(subseq))
#subseq = fa[chrom].antisense

#print(">{}/rc pyfastx\n{}\n".format(chrom, subseq))
