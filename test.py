import os
import time
import pyfastx

fasta = pyfastx.Fasta('test.fa.gz')

seq = fasta['seq10']

print(seq[1])


