import pyfastx

fa = pyfastx.Fasta('tests/data/test.fa.gz')

ids = fa.keys()

ids.sort(reverse=True)

ids.sort(key='name', reverse=True)