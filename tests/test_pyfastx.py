import os
import random
import pyfastx
import pyfaidx
import unittest

#os.chdir(os.path.dirname(__file__))

gzip_fasta = 'tests/data/test.fa.gz'
flat_fasta = 'tests/data/test.fa'

class FastaTest(unittest.TestCase):
	def setUp(self):
		self.fastx = pyfastx.Fasta(gzip_fasta)
		self.faidx = pyfaidx.Fasta(flat_fasta, sequence_always_upper=True)
		self.count = len(self.fastx)

	def tearDown(self):
		if os.path.exists('data/test.fa.gz.db'):
			os.remove('data/test.fa.gz.db')

		if os.path.exists('data/test.fa.fai'):
			os.remove('data/test.fa.fai')

	def get_random_index(self):
		return random.randint(0, self.count-1)

	def test_fasta(self):
		#seq counts
		self.assertEqual(len(self.fastx), len(self.faidx.keys()))

		#seq length
		expect = sum(len(s) for s in self.faidx)
		self.assertEqual(self.fastx.size, expect)

		expect = {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'N': 0}
		for s in self.faidx:
			expect['A'] += s.seq.count('A')
			expect['C'] += s.seq.count('C')
			expect['G'] += s.seq.count('G')
			expect['T'] += s.seq.count('T')
			expect['N'] += s.seq.count('N')

		#composition
		self.assertEqual(self.fastx.composition, expect)

		#GC content
		expect = (expect['G']+expect['C'])/sum(expect.values())

		self.assertEqual(self.fastx.gc_content, expect)


	def test_iter(self):
		for name, result in self.fastx:
			expect = self.faidx[name][:].seq
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
		self.assertEqual(expect[20:].seq, result[20:].seq)
		#self.assertEqual(expect[-10:-1].seq, result[-10:-1].seq)

	def test_seq_content(self):
		idx = self.get_random_index()
		result = self.fastx[idx]
		expect = self.faidx[idx]

		content = {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'N': 0}
		content['A'] += expect.seq.count('A')
		content['C'] += expect.seq.count('C')
		content['G'] += expect.seq.count('G')
		content['T'] += expect.seq.count('T')
		content['N'] += expect.seq.count('N')

		expect_gc = (content['G']+content['C'])/sum(content.values())

		self.assertEqual(result.composition, content)
		self.assertEqual(result.gc_content, expect_gc)

	def test_get_seq(self):
		idx = self.get_random_index()
		name = list(self.faidx.keys())[idx]
		l = len(self.fastx[idx])

		#test one interval
		a = int(l/2)
		interval = (random.randint(1, a), random.randint(a+1, l))

		expect = self.faidx.get_seq(name, interval[0], interval[1]).seq
		result = self.fastx.fetch(name, interval)

		self.assertEqual(expect, result)
