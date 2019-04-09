import os
import pyfastx
#os.remove('test.fa.gz.db')
fasta = pyfastx.fastx('test.fa.gz')
fasta.build_index()
print(fasta.get_sub_seq('seq1', 1, 2))
