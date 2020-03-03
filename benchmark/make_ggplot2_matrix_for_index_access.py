import sys

tools = ['pyfaidx', 'pyfasta', 'pysam', 'samtools', 'seqkit', 'pyfastx', 'pyfastx gzip']

species = {}
with open('species_name.txt') as fh:
	for line in fh:
		cols = line.strip().split('\t')
		species[cols[1]] = cols[0]

print("tool\tgenome\tsize\tcount\tmemory\ttime")

with open(sys.argv[1]) as fh:
	for line in fh:
		cols = line.strip().split()
		p = 0
		acc = species[cols[p]]
		p += 1
		size = str(round(int(cols[p])/1000000000, 5))
		p += 1
		count = cols[p]

		for i in range(7):
			p += 1
			memory = str(round(float(cols[p])/(1024*1024), 5))
			p += 1
			time = cols[p]

			print("{}\t{}\t{}\t{}\t{}\t{}".format(tools[i], acc, size, count, memory, time))

