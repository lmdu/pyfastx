import sys
import pyfastx

for read in pyfastx.Fastq(sys.argv[1]):
	print(read.name)
