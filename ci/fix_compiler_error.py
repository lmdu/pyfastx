import os
import sys

filename = os.path.join(sys.argv[1], r'Lib\distutils\cygwinccompiler.py')

with open(filename) as fp:
	if sys.version_info[0:2] > (3, 5):
		lines = fp.read().replace("            return ['msvcr100']", "            return ['msvcr100']\n        elif msc_ver == '1916':\n            return ['vcruntime140']")

	else:
		lines = fp.read().replace("            return ['msvcr100']", "            return ['msvcr100']\n        elif msc_ver == '1900':\n            return ['vcruntime140']")

with open(filename, "w") as fw:
	fw.write(lines)
