import os
import time
import pyfastx
#os.remove('test.fa.gz.db')
start = time.time()
fasta = pyfastx.Fasta('grcz10.fa')
print(time.time()-start)
#fasta.build_index()
#print(fasta.get_sub_seq('seq1', 1, 2))
