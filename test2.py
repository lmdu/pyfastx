import os
import time
import pyfaidx
#os.remove('test.fa.gz.db')
if os.path.exists('genome.fa.fai'):
	os.remove('genome.fa.fai')

start = time.time()
fasta = pyfaidx.Fasta('genome.fa')
print(time.time()-start)
