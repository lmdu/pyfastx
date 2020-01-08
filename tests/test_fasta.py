import os
import gzip
import random
import pyfastx
import pyfaidx
import unittest

gzip_fasta = 'tests/data/test.fa.gz'
flat_fasta = 'tests/data/test.fa'

class FastaTest(unittest.TestCase):
	def setUp(self):
		self.fastx = pyfastx.Fasta(gzip_fasta, build_index=False)
		self.fastx.build_index()
		self.fastx.rebuild_index()

		#reload index
		self.fastx = pyfastx.Fasta(gzip_fasta)

		self.fasta = pyfastx.Fasta(flat_fasta)

		self.faidx = pyfaidx.Fasta(flat_fasta, sequence_always_upper=True)
		
		self.count = len(self.fastx)

	def tearDown(self):
		if os.path.exists('{}.fxi'.format(gzip_fasta)):
			os.remove('{}.fxi'.format(gzip_fasta))

		if os.path.exists('{}.fxi'.format(flat_fasta)):
			os.remove('{}.fxi'.format(flat_fasta))

		if os.path.exists('{}.fai'.format(flat_fasta)):
			os.remove('{}.fai'.format(flat_fasta))

	def get_random_index(self):
		return random.randint(0, self.count-1)

	def test_module(self):
		# gzip check test
		self.assertEqual(pyfastx.gzip_check(gzip_fasta), self.fastx.is_gzip)

		# version test
		with open('src/version.h') as fh:
			version = fh.read().split()[2].strip('"')
			self.assertEqual(version, pyfastx.version())

		print(pyfastx.version(debug=True))

	def test_fasta(self):
		#fasta format
		self.assertEqual(self.fastx.type, 'DNA')

		#seq counts
		self.assertEqual(len(self.fastx), len(self.faidx.keys()))

		#seq length
		expect_size = sum(len(s) for s in self.faidx)
		self.assertEqual(self.fastx.size, expect_size)

		#test composition
		expect = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
		for s in self.faidx:
			expect['A'] += s[:].seq.count('A')
			expect['C'] += s[:].seq.count('C')
			expect['G'] += s[:].seq.count('G')
			expect['T'] += s[:].seq.count('T')
		self.assertEqual(self.fastx.composition, expect)

		#test GC content
		expect_gc = (expect['G']+expect['C'])/sum(expect.values())*100
		self.assertEqual(round(self.fastx.gc_content, 3), round(expect_gc, 3))

		#test GC skew
		expect_skew = (expect['G']-expect['C'])/(expect['G']+expect['C'])
		self.assertEqual(round(self.fastx.gc_skew, 3), round(expect_skew, 3))

		#test longest and shortest sequence
		longest = (None, 0)
		shortest = (None, expect_size)
		for seq in self.faidx:
			l = len(seq)
			if l > longest[1]:
				longest = (seq.name, l)

			if l < shortest[1]:
				shortest = (seq.name, l)

		long_seq = self.fastx.longest
		short_seq = self.fastx.shortest

		self.assertEqual(longest, (long_seq.name, len(long_seq)))
		self.assertEqual(shortest, (short_seq.name, len(short_seq)))

		#test contains
		idx = self.get_random_index()
		name = self.faidx[idx].name
		self.assertTrue(name in self.fastx)

	def test_iter_object(self):
		for seq in self.fastx:
			expect = self.faidx[seq.name][:].seq
			self.assertEqual(expect, seq.seq)

	def test_iter_tuple(self):
		fa = pyfastx.Fasta(gzip_fasta, build_index=False)
		
		for name, seq in fa:
			expect = str(self.faidx[name])
			self.assertEqual(expect, seq)

	def test_statistics(self):
		lens = sorted([len(seq) for seq in self.faidx], reverse=True)
		half = sum(lens)/2
		tlen = 0
		l50 = 0
		for n50 in lens:
			l50 += 1
			tlen += n50

			if tlen >= half:
				break

		self.assertEqual(self.fastx.nl(50), (n50, l50))

		#test mean length
		expect = round(sum(lens)/len(lens), 3)
		result = round(self.fastx.mean, 3)
		self.assertEqual(expect, result)

		#test median length
		lens = sorted(lens)
		expect = lens[105]

		result = self.fastx.median
		self.assertEqual(expect, result)

		#test count squence
		expect = 0
		for l in lens:
			if l >= 200:
				expect += 1
		result = self.fastx.count(200)
		self.assertEqual(expect, result)

	def test_key_identifier(self):
		fikeys = list(self.faidx.keys())
		fxkeys = list(self.fastx.keys())
		keyobj = self.fastx.keys()

		self.assertEqual(sorted(fikeys), sorted(fxkeys))

		#id counts
		self.assertEqual(len(fikeys), len(keyobj))

		idx = self.get_random_index()
		#get id from identifier class
		self.assertEqual(fikeys[idx], keyobj[idx])

		#negative index
		self.assertEqual(fikeys[len(fikeys)-idx], keyobj[len(fikeys)-idx])

		#check contains
		self.assertTrue(self.faidx[idx].name in keyobj)

	def test_keys_sort(self):
		#sort by id
		keys = self.fastx.keys()

		expect = [seq.name for seq in self.faidx]
		expect.reverse()
		result = [name for name in keys.sort('id', reverse=True)]
		self.assertEqual(expect, result)

		#sort by name
		expect = sorted([seq.name for seq in self.faidx])
		result = [name for name in keys.sort('name')]
		self.assertEqual(expect, result)

		#sort by length
		lens = [(seq.name, len(seq)) for seq in self.faidx]
		expect = [it[0] for it in sorted(lens, key=lambda x: x[1])]
		result = [name for name in keys.sort('length')]
		self.assertEqual(expect, result)

		keys.reset()

	def test_keys_filter(self):
		ids = self.fastx.keys()

		#test greater than
		expect = list(ids.filter(ids>700))
		result = [seq.name for seq in self.faidx if len(seq) > 700]
		self.assertEqual(expect, result)

		#test two compare
		expect = list(ids.filter(600<=ids<=700))
		result = [seq.name for seq in self.faidx if len(seq) >= 600 and len(seq) <= 700]
		self.assertEqual(expect, result)

		#test like compare
		expect = list(ids.filter(ids % 'JZ8226'))
		result = [seq.name for seq in self.faidx if seq.name.startswith('JZ8226')]
		self.assertEqual(expect, result)

		#test all compare
		expect = list(ids.filter(ids % 'JZ8226', ids>=300).sort('name', reverse=True))
		result = [seq.name for seq in self.faidx if seq.name.startswith('JZ8226') and len(seq) >= 300]
		self.assertEqual(sorted(result, reverse=True), expect)

		ids.reset()

	def test_seq_by_index(self):
		#test get seq by index
		idx = self.get_random_index()
		expect = self.faidx[idx][:]
		result = self.fastx[idx]

		self.assertEqual(expect.name, result.name)
		self.assertEqual(expect.seq, result.seq)

		#test subseq
		self.assertEqual(expect[0:10].seq, result[0:10].seq)

		#test negative index
		idx = (self.get_random_index() + 1) * -1
		expect = self.faidx[idx][:]
		result = self.fastx[idx]

		self.assertEqual(expect.name, result.name)
		self.assertEqual(expect.seq, result.seq)

	def test_seq_by_key(self):
		idx = self.get_random_index()
		key = list(self.faidx.keys())[idx]

		expect = self.faidx[key][:]
		result = self.fastx[key]

		self.assertEqual(expect.name, result.name)
		self.assertEqual(expect.seq, result.seq)

	def test_seq_reverse_complement(self):
		idx = self.get_random_index()
		expect = self.faidx[idx][:]
		result = self.fastx[idx]

		self.assertEqual(str(expect.reverse), result.reverse)
		self.assertEqual(str(expect.complement), result.complement)
		self.assertEqual(str(-expect), result.antisense)

	def test_seq_raw(self):
		idx = self.get_random_index()

		seq = self.fastx[idx]

		lines = []
		with open(flat_fasta) as fh:
			for line in fh:
				if line.startswith('>{}'.format(seq.name)):
					lines.append(line)
					continue

				if lines:
					if line[0] == '>':
						break

					lines.append(line)

		self.assertEqual(''.join(lines).replace('\n','\r\n'), seq.raw)

	def test_seq_slice(self):
		idx = self.get_random_index()
		expect = str(self.faidx[idx])
		result = self.fastx[idx]
		flatseq = self.fasta[idx]

		#test gzip subseq
		self.assertEqual(expect[5:10], result[5:10].seq)
		
		#test flat subseq
		self.assertEqual(expect[5:10], flatseq[5:10].seq)


		self.assertEqual(expect[20:], result[20:].seq)

		#test two level slice
		self.assertEqual(expect[10:100][:20], result[10:100][:20].seq)

		#test sequence index
		pos = random.randint(0, len(expect) - 1)
		self.assertEqual(expect[pos], result[pos])

		pos = random.randint(1, len(expect)) * -1
		self.assertEqual(expect[pos], result[pos])

		del flatseq

	def test_seq_content(self):
		idx = self.get_random_index()
		result = self.fastx[idx]
		expect = self.faidx[idx]

		content = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
		content['A'] += expect[:].seq.count('A')
		content['C'] += expect[:].seq.count('C')
		content['G'] += expect[:].seq.count('G')
		content['T'] += expect[:].seq.count('T')

		expect_gc = (content['G']+content['C'])/sum(content.values())*100

		self.assertEqual(result.composition, content)
		self.assertEqual(round(result.gc_content, 3), round(expect_gc, 3))

		#test gc skew
		expect_skew = (content['G']-content['C'])/(content['G']+content['C'])
		self.assertEqual(round(result.gc_skew, 3), round(expect_skew, 3))

	def test_seq_iter(self):
		idx = self.get_random_index()
		fai_seq = self.faidx[idx]
		fxi_seq = self.fastx[idx]
		fas_seq = self.fasta[idx]

		# test read seq line by line
		flatsq = [line for line in fas_seq]
		expect = [str(line) for line in fai_seq]
		result = [line for line in fxi_seq]

		self.assertEqual(expect, result)
		self.assertEqual(expect, flatsq)

		# test seq long name
		self.assertEqual(fai_seq.long_name.strip(), fxi_seq.description)

		# test seq str
		self.assertEqual(str(fai_seq), str(fxi_seq))

		# test seq contains
		s, e = sorted(random.sample(range(1, len(fai_seq)), 2))
		segment = str(fai_seq)[s-1:e]
		self.assertTrue(segment in fxi_seq)

		# test seq search
		expect = str(fai_seq).index(segment) + 1
		result = fxi_seq.search(segment)
		self.assertEqual(expect, result)

	def test_get_seq(self):
		idx = self.get_random_index()
		name = list(self.faidx.keys())[idx]
		l = len(self.fastx[idx])

		#test one interval
		a = int(l/2)
		interval = (random.randint(1, a), random.randint(a+1, l))

		expect = str(self.faidx[name])[interval[0]-1:interval[1]]
		result = self.fastx.fetch(name, interval)

		self.assertEqual(expect, result)

		#test multiple intervals
		intervals = []
		intervals.append((random.randint(1, int(a/2)), random.randint(int(a/2)+1, a)))
		intervals.append((random.randint(a+1, int((a+l)/2)), random.randint(int((a+l)/2)+1, l)))

		expect = "".join([str(self.faidx[name])[s-1:e] for s, e in intervals])
		result = self.fastx.fetch(name, intervals)

		self.assertEqual(expect, result)

if __name__ == '__main__':
	unittest.main()
