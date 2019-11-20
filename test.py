import pyfastx
fa=pyfastx.Fasta('tests/data/test.fa.gz')

ids=fa.keys()

for i in range(1, 200, 3):
	ids.filter(ids>i)
	print(ids)
