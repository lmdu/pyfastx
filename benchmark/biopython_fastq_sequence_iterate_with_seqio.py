import sys
from Bio import SeqIO

for seq in SeqIO.parse(sys.argv[1], 'fastq'):
	print(seq.seq)
