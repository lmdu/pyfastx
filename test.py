import pyfastx
fasta = pyfastx.fastx('test.fa.gz')

print(fasta.get_sub_seq('seq1', 1, 2))
