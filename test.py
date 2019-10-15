import pyfastx
fa = pyfastx.Fasta('test.fa')
fa = pyfastx.Fasta('test.fa', key_func=lambda x: x.split('|')[0])
fa.rebuild_index()
