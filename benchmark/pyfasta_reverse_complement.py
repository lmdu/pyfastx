import sys
import pyfasta

chrom, fafile = sys.argv[1:]

fa = pyfasta.Fasta(fafile, key_fn=lambda key: key.split()[0])

subseq = fa.sequence({'chr': chrom, 'start': 1, 'stop': len(fa[chrom]), 'strand': '-'})

print(">{}/rc pyfastx\n{}\n".format(chrom, subseq))
