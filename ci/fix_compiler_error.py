import os
import sys

filename = os.path.join(sys.argv[1], r'Lib\distutils\cygwinccompiler.py')

with open(filename) as fp:
	lines = f.read().replace("            return ['msvcr100']", "            return ['msvcr100']\n        elif msc_ver == '1900':\n            return ['vcruntime140']")

with open(filename, "w") as fw:
	fw.write(lines)
