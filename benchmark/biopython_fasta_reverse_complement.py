import sys
from Bio import SeqIO

longname, fafile = sys.argv[1:]

fa = SeqIO.index_db("{}.db".format(fafile))
s = fa[longname].reverse_complement()
print(str(s))
