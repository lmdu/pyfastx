import sys
import pysam
for read in pysam.FastxFile(sys.argv[1]):
	print(read.sequence)
