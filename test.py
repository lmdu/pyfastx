import os
import sys
import pyfastx

fqfile = sys.argv[1]
if os.path.exists('{}.fxi'.format(fqfile)):
    os.remove('{}.fxi'.format(fqfile))
fq=pyfastx.Fastq(fqfile)
