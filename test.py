import pyfastx

fa = pyfastx.Fasta('hg19.fa')

start0 = 60757997
end0 = 60758118
for i in range(0,4):
	print('start='+str(start0+i)+', end='+str(end0+i))
	print(str(fa['chr20'][start0+i:end0+i]))
