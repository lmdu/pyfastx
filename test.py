import os
import sys
import pyfastx

fqfile = 'tests/data/test.fa.gz'
if os.path.exists('{}.fxi'.format(fqfile)):
    os.remove('{}.fxi'.format(fqfile))
fa=pyfastx.Fasta(fqfile)

print(len(fa))

s = fa[0]

print(s.name)
