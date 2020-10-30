import sys
from Bio import SeqIO
SeqIO.index_db("{}.db".format(sys.argv[1]), sys.argv[1], 'fasta')
