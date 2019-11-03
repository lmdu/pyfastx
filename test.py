import os
import time
import pyfaidx
import pyfastx

fa_file = '/mnt/d/research/tandem/data/t.fa'

if os.path.exists('%s.fai' % fa_file):
	os.remove('%s.fai' % fa_file)

start = time.time()
fa = pyfaidx.Fasta(fa_file)
print(time.time() - start)

if os.path.exists('%s.fxi' % fa_file):
	os.remove('%s.fxi' % fa_file)

start = time.time()
fx = pyfastx.Fasta(fa_file)
print(time.time() - start)