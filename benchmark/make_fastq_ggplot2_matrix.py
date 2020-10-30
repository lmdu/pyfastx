import sys
import csv

title = ["file", "fsize", "gsize", "count", "tool", "memory", "time"]

def kb2mb(v):
	return str(float(v)/1024)

def bp2gb(v):
	return str(float(v)/1000000000)

def kb2gb(v):
	return str(float(v)/1048576)

def bs2gb(v):
	return str(float(v)/1073741824)

with open(sys.argv[1]) as fh:
	reader = csv.reader(fh, delimiter='\t')
	header = next(reader)

	if not header[7]:
		title.append("index")

	print("\t".join(title))
	
	column = len(header)
	for row in reader:
		fn = row[0]
		fs = bs2gb(row[3])
		gs = bs2gb(row[4])
		c = row[2]

		i = 5
		while (i < column):
			items = [fn, fs, gs, c]
			items.append(header[i])
			items.append(kb2gb(row[i]))
			i += 1
			items.append(row[i])
			i += 1

			if not header[7]:
				items.append(kb2gb(row[i]))
				i += 1

			print("\t".join(items))
