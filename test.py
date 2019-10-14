import pyfastx
fa = pyfastx.Fasta('test.fa', key_func=lambda x: x.split('|')[1])