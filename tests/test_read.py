import os
import gzip
import random
import pyfastx
import unittest

join = os.path.join
data_dir = join(os.path.dirname(os.path.abspath(__file__)), 'data')

gzip_fastq = join(data_dir, 'test.fq.gz')
flat_fastq = join(data_dir, 'test.fq')

class ReadTest(unittest.TestCase):
	def setUp(self):

		self.fastq = pyfastx.Fastq(gzip_fastq)

		#reload index
		self.fastq = pyfastx.Fastq(gzip_fastq)

		#flat fastq
		self.flatq = pyfastx.Fastq(flat_fastq)

		self.reads = {}
		self.bases = {'A': 0, 'T': 0, 'G': 0, 'C':0, 'N':0}
		i = 0
		c = -1
		with open(flat_fastq) as fh:
			for line in fh:
				i += 1

				if i % 4 == 1:
					c += 1
					self.reads[c] = [line[1:].strip().split()[0], 0, 0]

				elif i % 4 == 2:
					self.reads[c][1] = line.strip()

					self.bases['A'] += line.count('A')
					self.bases['T'] += line.count('T')
					self.bases['G'] += line.count('G')
					self.bases['C'] += line.count('C')
					self.bases['N'] += line.count('N')

				elif i % 4 == 0:
					self.reads[c][2] = line.strip()

	def tearDown(self):
		del self.fastq
		del self.flatq

		if os.path.exists('{}.fxi'.format(gzip_fastq)):
			os.remove('{}.fxi'.format(gzip_fastq))

		if os.path.exists('{}.fxi'.format(flat_fastq)):
			os.remove('{}.fxi'.format(flat_fastq))

	def get_random_read(self):
		return random.randint(0, len(self.fastq)-1)

	def test_repr(self):
		idx = self.get_random_read()
		read = self.flatq[idx]
		read1 = self.reads[idx]

		self.assertEqual(repr(read), "<Read> {} with length of {}".format(read1[0], len(read1[1])))

	def test_read(self):
		idx = self.get_random_read()
		result = self.fastq[idx]
		expect = self.reads[idx]

		del result
		result = self.fastq[idx]

		read0 = self.flatq[idx]

		# test length
		self.assertEqual(len(result), len(expect[1]))

		# test name
		self.assertEqual(result.name, expect[0])

		# test str
		self.assertEqual(str(result), expect[1])

		# test seq
		self.assertEqual(result.seq, expect[1])
		self.assertEqual(read0.seq, expect[1])

		# test quality
		self.assertEqual(result.qual, expect[2])
		self.assertEqual(read0.qual, expect[2])

		# test quality integer
		self.assertEqual(result.quali, [ord(b)-33 for b in expect[2]])

		result = self.fastq[expect[0]]

		# test subscript
		self.assertEqual(result.seq, expect[1])

		# test contain
		self.assertTrue(result.name in self.fastq)

	def test_read_seq(self):
		#test reverse
		idx = self.get_random_read()
		result = self.fastq[idx]
		expect = self.reads[idx][1]
		self.assertEqual(result.reverse, expect[::-1])

		#test complement
		complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N':'N'}
		expect1 = ''.join(complement.get(b) for b in expect)
		self.assertEqual(result.complement, expect1)

		#test antisense
		expect2 = ''.join(complement.get(b) for b in reversed(expect))
		self.assertEqual(result.antisense, expect2)

	def test_read_description(self):
		idx = self.get_random_read()
		read = self.fastq[idx]

		i = -1
		with open(flat_fastq) as fh:
			for line in fh:
				if line[0] == '@':
					i += 1

					if i == idx:
						break

		self.assertEqual(line.strip(), read.description)

		read = self.flatq[idx]
		self.assertEqual(line.strip(), read.description)

	def test_read_raw(self):
		idx = self.get_random_read()
		read = self.fastq[idx]

		lines = []
		with gzip.open(gzip_fastq, 'rb') as fh:
			for line in fh:
				line = line.decode()
				if line.startswith('@{}'.format(read.name)):
					lines.append(line)
					continue

				if lines:
					if line[0] == '@':
						break

					lines.append(line)

		self.assertEqual(''.join(lines), read.raw)

		idx = self.get_random_read()
		read = self.flatq[idx]
		lines = []
		with open(flat_fastq, 'rb') as fh:
			for line in fh:
				line = line.decode()
				if line.startswith('@{}'.format(read.name)):
					lines.append(line)
					continue

				if lines:
					if line[0] == '@':
						break

					lines.append(line)

		self.assertEqual(''.join(lines), read.raw)

if __name__ == '__main__':
	unittest.main()
