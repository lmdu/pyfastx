import os
import pyfaidx
import pyfastx
import unittest

join = os.path.join
data_dir = join(os.path.dirname(os.path.abspath(__file__)), 'data')

gzip_fasta = join(data_dir, 'test.fa.gz')
gzip_fastq = join(data_dir, 'test.fq.gz')
flat_fasta = join(data_dir, 'test.fa')
flat_fastq = join(data_dir, 'test.fq')

class FastxTest(unittest.TestCase):
	def setUp(self):
		self.faidx = pyfaidx.Fasta(flat_fasta)

	def tearDown(self):
		del self.faidx

		if os.path.exists('{}.fai'.format(flat_fasta)):
			os.remove('{}.fai'.format(flat_fasta))

	def test_fasta_iter(self):
		for name, seq, comment in pyfastx.Fastx(gzip_fasta, comment=True):
			s = self.faidx[name]
			self.assertEqual(str(seq), seq)
			self.assertEqual(' '.join(s.long_name.split()[1:]), comment)

	def test_fasta_upper(self):
		for name, seq in pyfastx.Fastx(flat_fasta, uppercase=True):
			self.assertEqual(str(self.faidx[name]), seq)

	def test_fastq_iter(self):
		reads = {}
		with open(flat_fastq) as fh:
			for line in fh:
				line = line.strip()

				if line[0] == '@':
					name = line.split()[0][1:]
					comment = ' '.join(line.split()[1:])
					reads[name] = [comment]
				elif line == '+':
					continue
				else:
					reads[name].append(line)

		for name, seq, qual, comment in pyfastx.Fastx(gzip_fastq, "fastq", comment=True):
			r = reads[name]

			self.assertEqual(r[0], comment)
			self.assertEqual(r[1], seq)
			self.assertEqual(r[2], qual)

	def test_fastx_repr(self):
		fa = pyfastx.Fastx(gzip_fasta, "fasta")
		self.assertEqual(repr(fa), "<Fastx> fasta {}".format(gzip_fasta))

		fq = pyfastx.Fastx(gzip_fastq, "fastq")
		self.assertEqual(repr(fq), "<Fastx> fastq {}".format(gzip_fastq))

	def test_exception(self):
		with self.assertRaises(FileExistsError):
			_ = pyfastx.Fastx('test_file')

		with self.assertRaises(RuntimeError):
			_ = pyfastx.Fastx(gzip_fasta, format="fastx")
