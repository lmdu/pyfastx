import os
import random
import pyfastx
import pyfaidx
import unittest

join = os.path.join
data_dir = join(os.path.dirname(os.path.abspath(__file__)), 'data')

flat_fasta = join(data_dir, 'test.fa')

class SequenceErrorTest(unittest.TestCase):
	def setUp(self):
		self.fastx = pyfastx.Fasta(flat_fasta)

		self.faidx = pyfaidx.Fasta(flat_fasta, sequence_always_upper=True)

		self.count = len(self.fastx)

	def tearDown(self):
		del self.fastx
		del self.faidx

		if os.path.exists('{}.fxi'.format(flat_fasta)):
			os.remove('{}.fxi'.format(flat_fasta))

		if os.path.exists('{}.fai'.format(flat_fasta)):
			os.remove('{}.fai'.format(flat_fasta))

	def get_random_index(self):
		return random.randint(0, self.count-1)

	def get_random_key(self):
		idx = self.get_random_index()
		return list(self.faidx.keys())[idx]

	def test_seq_by_index(self):
		# get valid index
		idx = self.get_random_index()

		expect = self.faidx[idx][:]
		result = self.fastx[idx]

		self.assertEqual(expect.name, result.name)
		self.assertEqual(expect.seq, result.seq)

		# get invalid index
		with self.assertRaises(IndexError):
		    self.fastx[self.count + 100]

		# get valid index again
		expect = self.faidx[idx][:]
		result = self.fastx[idx]

		self.assertEqual(expect.name, result.name)
		self.assertEqual(expect.seq, result.seq)

	def test_seq_by_key(self):
		# get valid key
		key = self.get_random_key()

		expect = self.faidx[key][:]
		result = self.fastx[key]

		self.assertEqual(expect.name, result.name)
		self.assertEqual(expect.seq, result.seq)

		# get invalid key
		with self.assertRaises(KeyError):
		    self.fastx['__BOGUS_KEY_FOR_TESTING__']

		# get valid key again
		expect = self.faidx[key][:]
		result = self.fastx[key]

		self.assertEqual(expect.name, result.name)
		self.assertEqual(expect.seq, result.seq)


if __name__ == '__main__':
	unittest.main()
