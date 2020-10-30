import sys
import pyfastx

fq = pyfastx.Fastq(sys.argv[1])
print(fq.size, len(fq))
