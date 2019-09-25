import os
import random
import pyfastx
import pyfaidx
import unittest
import statistics

#os.chdir(os.path.dirname(__file__))

gzip_fasta = 'tests/data/test.fa.gz'
flat_fasta = 'tests/data/test.fa'

class FastaTest(unittest.TestCase):
	def setUp(self):
		self.fastx = pyfastx.Fasta(gzip_fasta)
		self.faidx = pyfaidx.Fasta(flat_fasta, sequence_always_upper=True)
		self.count = len(self.fastx)

	def tearDown(self):
		if os.path.exists('tests/data/test.fa.gz.fxi'):
			os.remove('tests/data/test.fa.gz.fxi')

		if os.path.exists('tests/data/test.fa.fai'):
			os.remove('tests/data/test.fa.fai')

	def get_random_index(self):
		return random.randint(0, self.count-1)

	def test_fasta(self):
		#seq counts
		self.assertEqual(len(self.fastx), len(self.faidx.keys()))

		#seq length
		expect_size = sum(len(s) for s in self.faidx)
		self.assertEqual(self.fastx.size, expect_size)

		#test composition
		expect = {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'N': 0}
		for s in self.faidx:
			expect['A'] += s[:].seq.count('A')
			expect['C'] += s[:].seq.count('C')
			expect['G'] += s[:].seq.count('G')
			expect['T'] += s[:].seq.count('T')
			expect['N'] += s[:].seq.count('N')
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

		self.assertEqual(longest, self.fastx.longest)
		self.assertEqual(shortest, self.fastx.shortest)

		#test contains
		idx = self.get_random_index()
		name = self.faidx[0].name
		self.assertTrue(name in self.fastx)

	def test_iter(self):
		for name, result in self.fastx:
			expect = self.faidx[name][:].seq
			self.assertEqual(expect, result)

	def test_stat(self):
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
		expect = round(statistics.mean(lens), 3)
		result = round(self.fastx.mean, 3)
		self.assertEqual(expect, result)

		#test median length
		expect = statistics.median(lens)
		result = self.fastx.median
		self.assertEqual(expect, result)

		#test count squence
		expect = 0
		for l in lens:
			if l >= 200:
				expect += 1
		result = self.fastx.count(200)
		self.assertEqual(expect, result)

	def test_keys(self):
		expect = list(self.faidx.keys())
		result = list(self.fastx.keys())

		self.assertEqual(expect, result)

		#id counts
		ids = self.fastx.keys()
		self.assertEqual(len(ids), len(expect))

		#get id from identifier class
		self.assertEqual(ids[0], expect[0])
		self.assertEqual(ids[-1], expect[-1])

		#check contains
		idx = self.get_random_index()
		name = self.faidx[idx].name

		self.assertTrue(name in ids)

	def test_seq_by_index(self):
		#test get seq by index
		idx = self.get_random_index()
		expect = self.faidx[idx][:]
		result = self.fastx[idx]

		self.assertEqual(expect.name, result.name)
		self.assertEqual(expect.seq, result.seq)

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

		self.assertEqual(expect.reverse.seq, result.reverse)
		self.assertEqual(expect.complement.seq, result.complement)
		expect = -expect
		self.assertEqual(expect.seq, result.antisense)

	def test_seq_slice(self):
		idx = self.get_random_index()
		expect = self.faidx[idx]
		result = self.fastx[idx]

		self.assertEqual(expect[5:10].seq, result[5:10].seq)
		#self.assertEqual(expect[20:].seq, result[20:].seq)
		expect = expect[20:].seq
		result = result[20:].seq
		if len(expect) > len(result):
			expect = expect[0:-1]
		
		self.assertEqual(expect, result)

	def test_seq_content(self):
		idx = self.get_random_index()
		result = self.fastx[idx]
		expect = self.faidx[idx]

		content = {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'N': 0}
		content['A'] += expect[:].seq.count('A')
		content['C'] += expect[:].seq.count('C')
		content['G'] += expect[:].seq.count('G')
		content['T'] += expect[:].seq.count('T')
		content['N'] += expect[:].seq.count('N')

		expect_gc = (content['G']+content['C'])/sum(content.values())*100

		self.assertEqual(result.composition, content)
		self.assertEqual(round(result.gc_content, 3), round(expect_gc, 3))

	def test_get_seq(self):
		idx = self.get_random_index()
		name = list(self.faidx.keys())[idx]
		l = len(self.fastx[idx])

		#test one interval
		a = int(l/2)
		interval = (random.randint(1, a), random.randint(a+1, l))

		#expect = self.faidx.get_seq(name, interval[0], interval[1]).seq
		#expect = self.faidx[name][interval[0]-1:interval[1]].seq
		expect = str(self.faidx[name])[interval[0]-1:interval[1]]
		result = self.fastx.fetch(name, interval)

		if len(expect) > len(result):
			expect = expect[0:-1]

		self.assertEqual(expect, result)
