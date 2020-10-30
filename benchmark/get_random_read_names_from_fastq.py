import sys
import random
import pyfastx

random.seed(sys.argv[1])

for fqfile in sys.argv[2:]:
	fq = pyfastx.Fastq(fqfile)
	samples = set(random.sample(range(len(fq)), 10000))
	with open("{}.list".format(fqfile), 'w') as fw:
		for r in fq:
			if (r.id - 1) in samples: 
				fw.write("{}\n".format(r.name))

	print(fqfile)
	
