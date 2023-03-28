import os
import random
import pyfastx
import unittest

join = os.path.join
data_dir = join(os.path.dirname(os.path.abspath(__file__)), 'data')

gzip_fastq = join(data_dir, 'test.fq.gz')
flat_fastq = join(data_dir, 'test.fq')

class FastqTest(unittest.TestCase):
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

	def test_build(self):
		del self.fastq

		if os.path.exists('{}.fxi'.format(gzip_fastq)):
			os.remove('{}.fxi'.format(gzip_fastq))

		fq = pyfastx.Fastq(gzip_fastq, build_index=False)
		fq.build_index()

		self.fastq = pyfastx.Fastq(gzip_fastq)

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
		self.assertEqual(['Sanger Phred+33', 'Illumina 1.8+ Phred+33', 'PacBio HiFi Phred+33'], self.fastq.encoding_type)

		# test phred
		self.assertEqual(self.fastq.phred, 33)

		#test length
		self.assertEqual(self.fastq.avglen, 150)
		self.assertEqual(self.fastq.minlen, 150)
		self.assertEqual(self.fastq.maxlen, 150)

		#test quality
		self.assertEqual(self.fastq.maxqual, 70)
		self.assertEqual(self.fastq.minqual, 35)

	def test_platform(self):
		#test unknown
		fqfile = 'test1.fq'
		with open(fqfile, 'w') as fw:
			fw.write("@read1\n")
			fw.write("AAAAAAAAAA\n")
			fw.write("+\n")
			fw.write("          \n")
		fq = pyfastx.Fastq(fqfile)
		self.assertEqual(fq.encoding_type, ['Unknown'])
		del fq
		os.remove("{}.fxi".format(fqfile))

		with open(fqfile, 'w') as fw:
			fw.write("@read1\n")
			qs = []
			for i in range(59, 105):
				qs.append(chr(i))
			fw.write("{}\n".format('A'*len(qs)))
			fw.write("+\n")
			fw.write("{}\n".format(''.join(qs)))
		fq = pyfastx.Fastq(fqfile)
		self.assertIn("Solexa Solexa+64", fq.encoding_type)
		del fq
		os.remove("{}.fxi".format(fqfile))

		with open(fqfile, 'w') as fw:
			fw.write("@read1\n")
			qs = []
			for i in range(64, 105):
				qs.append(chr(i))
			fw.write("{}\n".format('A'*len(qs)))
			fw.write("+\n")
			fw.write("{}\n".format(''.join(qs)))
		fq = pyfastx.Fastq(fqfile)
		self.assertIn("Illumina 1.3+ Phred+64", fq.encoding_type)
		del fq
		os.remove("{}.fxi".format(fqfile))

		with open(fqfile, 'w') as fw:
			fw.write("@read1\n")
			qs = []
			for i in range(66, 105):
				qs.append(chr(i))
			fw.write("{}\n".format('A'*len(qs)))
			fw.write("+\n")
			fw.write("{}\n".format(''.join(qs)))
		fq = pyfastx.Fastq(fqfile)
		self.assertIn("Illumina 1.5+ Phred+64", fq.encoding_type)
		del fq
		os.remove("{}.fxi".format(fqfile))
		os.remove(fqfile)

	def test_negative(self):
		read = self.fastq[-1]
		self.assertEqual(read.name, self.reads[len(self.reads)-1][0])

	def test_iter_object(self):
		# test read iter
		i = -1
		for read in self.fastq:
			i += 1
			self.assertEqual(read.name, self.reads[i][0])
			self.assertEqual(read.seq, self.reads[i][1])
			self.assertEqual(read.qual, self.reads[i][2])

		#test reference of read made from loop
		for read in self.fastq:
			break

		expect = self.reads[0][1]
		self.assertEqual(expect, read.seq)
		self.assertEqual(expect, read.seq)

	def test_iter_tuple(self):
		i = -1
		for name, seq, qual in pyfastx.Fastq(flat_fastq, build_index=False):
			i += 1
			self.assertEqual(name, self.reads[i][0])
			self.assertEqual(seq, self.reads[i][1])
			self.assertEqual(qual, self.reads[i][2])

	def test_read_len(self):
		lens = [len(it[1]) for it in self.reads.values()]

		self.assertEqual(self.fastq.minlen, min(lens))
		self.assertEqual(self.fastq.maxlen, min(lens))

	def test_repr(self):
		self.assertEqual(repr(self.fastq), "<Fastq> {} contains {} reads".format(gzip_fastq, len(self.reads)))

	def test_full_name(self):
		fq = pyfastx.Fastq(flat_fastq, build_index=False, full_name=True)

		for name, _, _ in fq:
			self.assertTrue(name, self.fastq[name.split()[0]].description)

	def test_exception(self):
		with self.assertRaises(FileExistsError):
			_ = pyfastx.Fastq('a_fastq_file_not_exists')

		with self.assertRaises(IndexError):
			_ = self.fastq[len(self.fastq)]

		with self.assertRaises(KeyError):
			_ = self.fastq[int]

		with self.assertRaises(KeyError):
			_ = self.fastq['abc']

if __name__ == '__main__':
	unittest.main()
