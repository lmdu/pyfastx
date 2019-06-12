import os
import time
import pyfastx
#os.remove('test.fa.gz.db')
if os.path.exists('genome.fa.db'):
	os.remove('genome.fa.db')

start = time.time()
fasta = pyfastx.Fasta('genome.fa.gz')
print(time.time()-start)
#fasta.build_index()
#print(fasta.get_sub_seq('seq1', 1, 2))
