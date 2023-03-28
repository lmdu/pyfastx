import os
import random
import pyfastx
import unittest

join = os.path.join
data_dir = join(os.path.dirname(os.path.abspath(__file__)), 'data')

gzip_fastq = join(data_dir, 'test.fq.gz')
flat_fastq = join(data_dir, 'test.fq')

class FastxTest(unittest.TestCase):
	def setUp(self):
		self.fastq = pyfastx.Fastq(gzip_fastq)

		with open(flat_fastq) as fh:
			self.keys = [line.split()[0][1:] for line in fh if line[0] == '@']
		self.count = len(self.keys)

	def tearDown(self):
		del self.fastq

		if os.path.exists('{}.fxi'.format(gzip_fastq)):
			os.remove('{}.fxi'.format(gzip_fastq))

	def get_random_index(self):
		return random.randint(0, self.count-1)

	def test_fastq_key(self):
		#test for repr
		keys = self.fastq.keys()
		self.assertEqual(repr(keys), "<FastqKeys> contains {} keys".format(self.count))

		#test for lenth
		self.assertEqual(self.count, len(keys))

		#test for iter
		for i, k in enumerate(keys):
			self.assertEqual(self.keys[i], k)

		#test random access
		idx = self.get_random_index()
		self.assertEqual(self.keys[idx], keys[idx])

		#test contain
		idx = self.get_random_index()
		self.assertTrue(self.keys[idx] in keys)

		#test negative index
		neg = (self.count - idx) * -1
		self.assertEqual(self.keys[neg], keys[neg])

	def test_exception(self):
		with self.assertRaises(IndexError):
			_ = self.keys[self.count*2]
