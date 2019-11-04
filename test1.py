import pyfastx
fa = pyfastx.Fasta('tests/data/test.fa.gz')
print(fa.gc_content)
print(fa.composition)
