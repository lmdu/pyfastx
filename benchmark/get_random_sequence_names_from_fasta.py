"""randomly get numbers of sequence name from fasta file
Usage:
	python3 get_random_sequence_names.py seed prop fastas...
@seed, random seed
@prop, proportion of sequence to extract, 0-1
"""

import sys
import math
import random
import pyfastx

random.seed(sys.argv[1])
prop = float(sys.argv[2])

for fafile in sys.argv[3:]:
	fa = pyfastx.Fasta(fafile)
	num = math.ceil(len(fa)*prop)
	ids = fa.keys()

	samples = random.sample(range(len(fa)), num)

	with open("{}.list".format(fafile), 'w') as fw:
		for i in samples:
			fw.write("{}\n".format(ids[i]))
