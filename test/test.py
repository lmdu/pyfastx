import os
import random
import pyfastx
import pyfaidx
import unittest

os.chdir(os.path.dirname(__file__))

gzip_fasta = 'data/test.fa.gz'
flat_fasta = 'data/test.fa'

class FastaTest(unittest.TestCase):
	def setUp(self):
		self.fastx = pyfastx.Fasta(gzip_fasta)
		self.faidx = pyfaidx.Fasta(flat_fasta, sequence_always_upper=True)
		self.count = len(self.faidx)

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

	def test_keys(self):
		expect = list(self.faidx.keys())
		result = list(self.fastx.keys())

		self.assertEqual(expect, result)

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
		self.assertEqual(expect.name, result.seq)

	def test_seq_reverse_complement(self):
		idx = self.get_random_index()
		expect = self.faidx[idx][:]
		result = self.fastx[idx]

		self.assertEqual(expect.reverse.seq, result.reverse)
		self.assertEqual(expect.complement.seq, result.complement)
		self.assertEqual(-expect, result.antisense)

	def test_seq_slice(self):
		idx = self.get_random_index()
		expect = self.faidx[idx]
		result = self.fastx[idx]

		self.assertEqual(expect[5:10].seq, result[5:10].seq)
		self.assertEqual(expect[20:].seq, result[20:].seq)
		self.assertEqual(expect[-10:-1].seq, result[-10:-1].seq)

	def test_get_seq(self):
		idx = self.get_random_index()
		expect = self.faidx[idx]
		result = self.fastx[idx]

		










