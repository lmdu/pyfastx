import os
import random
import pyfastx
import pyfaidx
import unittest

join = os.path.join
data_dir = join(os.path.dirname(os.path.abspath(__file__)), 'data')

gzip_fasta = join(data_dir, 'test.fa.gz')
flat_fasta = join(data_dir, 'test.fa')

class IdentifierTest(unittest.TestCase):
	def setUp(self):
		self.fastx = pyfastx.Fasta(gzip_fasta)

		self.fasta = pyfastx.Fasta(flat_fasta)

		self.faidx = pyfaidx.Fasta(flat_fasta, sequence_always_upper=True)

		self.count = len(self.fastx)

	def tearDown(self):
		del self.fastx
		del self.fasta
		del self.faidx

		if os.path.exists('{}.fxi'.format(gzip_fasta)):
			os.remove('{}.fxi'.format(gzip_fasta))

		if os.path.exists('{}.fxi'.format(flat_fasta)):
			os.remove('{}.fxi'.format(flat_fasta))

		if os.path.exists('{}.fai'.format(flat_fasta)):
			os.remove('{}.fai'.format(flat_fasta))

	def get_random_index(self):
		return random.randint(0, self.count-1)

	def test_key_repr(self):
		ids = self.fastx.keys()
		self.assertEqual(repr(ids), "<FastaKeys> contains {} keys".format(len(ids)))

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
		self.assertEqual(fikeys[idx-len(fikeys)], keyobj[idx-len(fikeys)])
		self.assertEqual(fikeys[-1], keyobj[-1])

		#check contains
		self.assertTrue(self.faidx[idx].name in keyobj)

	def test_key_slice(self):
		fikeys = list(self.faidx.keys())
		fxkeys = self.fastx.keys()
		idx = self.get_random_index()

		#test positive
		self.assertEqual(fikeys[idx:idx+10], fxkeys[idx:idx+10])

		#test negative
		idx = self.get_random_index() - self.count
		self.assertEqual(fikeys[idx:idx+10], fxkeys[idx:idx+10])

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
		result = sorted(result, reverse=True)
		self.assertEqual(result, expect)

		#get filtered item
		self.assertEqual(result[0], ids[0])
		self.assertEqual(result[-1], ids[-1])

		ids.reset()
		ids.filter(ids % 'JZ8226', ids >= 200)

	def test_id_exception(self):
		keyobj = self.fastx.keys()
		with self.assertRaises(IndexError):
			_ = keyobj[len(keyobj)]

		with self.assertRaises(ValueError):
			_ = keyobj % list

		with self.assertRaises(ValueError):
			_ = keyobj.filter()

		with self.assertRaises(ValueError):
			_ = keyobj.sort('sort')

		with self.assertRaises(ValueError):
			_ = keyobj > list

		with self.assertRaises(TypeError):
			_ = keyobj[list]

if __name__ == '__main__':
	unittest.main()
