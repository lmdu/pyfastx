import os
import gzip
import random
import pyfastx
import pyfaidx
import unittest

gzip_fastq = 'tests/data/test.fq.gz'
flat_fastq = 'tests/data/test.fq'

class FastaTest(unittest.TestCase):
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
		if os.path.exists('{}.fxi'.format(gzip_fastq)):
			os.remove('{}.fxi'.format(gzip_fastq))

		if os.path.exists('{}.fxi'.format(flat_fastq)):
			os.remove('{}.fxi'.format(flat_fastq))

	def get_random_read(self):
		return random.randint(0, len(self.fastq)-1)

	def test_build(self):
		self.fastx = None

		if os.path.exists('{}.fxi'.format(gzip_fastq)):
			os.remove('{}.fxi'.format(gzip_fastq))

		fq = pyfastx.Fastq(gzip_fastq, build_index=False)
		fq.build_index()
		fq = pyfastx.Fastq(gzip_fastq)

	def test_fastq(self):
		# test gzip format
		self.assertEqual(pyfastx.gzip_check(gzip_fastq), self.fastq.is_gzip)

		# test seq length
		self.assertEqual(self.fastq.size, sum(self.bases.values()))

		# test length
		self.assertEqual(len(self.reads), len(self.fastq))

		# test gc content
		result = round(self.fastq.gc_content, 2)
		
		expect = round((self.bases['G']+self.bases['C'])/(sum(self.bases.values()))*100, 2)
		self.assertEqual(expect, result)

		# test composition
		self.assertEqual(self.fastq.composition, self.bases)

		# test encoding type
		self.assertEqual(['Sanger Phred+33', 'Illumina 1.8+ Phred+33'], self.fastq.encoding_type)

		# test phred
		self.assertEqual(self.fastq.phred, 33)

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

	def test_iter_object(self):
		# test read iter
		i = -1
		for read in self.fastq:
			i += 1
			self.assertEqual(read.name, self.reads[i][0])
			self.assertEqual(read.seq, self.reads[i][1])
			self.assertEqual(read.qual, self.reads[i][2])

	def test_iter_tuple(self):
		i = -1
		for name, seq, qual in pyfastx.Fastq(flat_fastq, build_index=False):
			i += 1
			self.assertEqual(name, self.reads[i][0])
			self.assertEqual(seq, self.reads[i][1])
			self.assertEqual(qual, self.reads[i][2])

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


	def test_read_raw(self):
		idx = self.get_random_read()
		read = self.fastq[idx]

		i = -1
		lines = []
		with gzip.open(gzip_fastq, 'rt') as fh:
			for line in fh:
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
