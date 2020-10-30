import sys
import math
import random
import pyfastx

random.seed(sys.argv[1])

for fafile in sys.argv[2:]:
	fa = pyfastx.Fasta(fafile)
	ids = fa.keys()
	names = list(ids.filter(ids>1000))

	list_fh = open("{}.interval".format(fafile), 'w')
	region_fh = open("{}.region".format(fafile), 'w')

	if len(names) >= 1000:
		lists = random.sample(names, 1000)
	else:
		lists = random.sample(names * math.ceil(1000/len(names)), 1000)

	for name in lists:
		pos = random.randint(1, len(fa[name])-1000)
		list_fh.write("{}\t{}\t{}\n".format(name, pos, pos+1000))
		region_fh.write("{}:{}-{}\n".format(name, pos, pos+1000))

	list_fh.close()
	region_fh.close()
