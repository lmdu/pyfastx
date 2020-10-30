import sys
import csv

mapping = {
	'GCA_002915635.1': "Axolotl",
	'GCA_001447015.2': "Sugar pine",
	'GCA_000411955.5': "White spruce",
	'GCA_900519105.1': "Bread wheat",
	'GCA_900184675.1': "Wild wheat",
	'GCA_002284835.2': "Bullfrog",
	'GCA_000516895.1': "Migratory locust",
	'GCA_003013575.1': "Pea",
	'GCF_000001405.38': "Human",
	'GCF_000005005.2': "Maize",
	'GCF_000002315.5': "Chicken",
	'GCF_000003195.3': "Sorghum",
	'GCF_000331145.1': "Chickpea",
	'GCF_000001735.4': "Thale cress",
	'GCF_000146045.2': "Yeast",
	'SRR6649487_10': '10%',
	'SRR6649487_50': '50%',
	'SRR6649487_1': '100%'
}

title = ["genome", "size", "count", "tool", "memory", "time"]

def kb2mb(v):
	return str(float(v)/1024)

def bp2gb(v):
	return str(float(v)/1000000000)

def kb2gb(v):
	return str(float(v)/1048576)

with open(sys.argv[1]) as fh:
	reader = csv.reader(fh, delimiter='\t')
	header = next(reader)

	if not header[5]:
		title.append("index")

	print("\t".join(title))
	
	column = len(header)
	for row in reader:
		g = mapping.get(row[0], row[0])
		s = bp2gb(row[1])
		c = row[2]

		i = 3
		while (i < column):
			items = [g, s, c]
			items.append(header[i])
			items.append(kb2gb(row[i]))
			i += 1
			items.append(row[i])
			i += 1

			if not header[5]:
				items.append(kb2mb(row[i]))
				i += 1

			print("\t".join(items))
