import os
import gzip
import shutil
import pyfastx
import tempfile
import unittest

join = os.path.join

COMMENTED = ";comment line\n;~key: value\n;\n>chr1 first\nACGTACGT\nACGT\n>chr2\nTTTTGGGG\n"


class CommentTest(unittest.TestCase):
	def setUp(self):
		self.tmp_dir = tempfile.mkdtemp()

	def tearDown(self):
		shutil.rmtree(self.tmp_dir, ignore_errors=True)

	def write(self, name, text, compress=False):
		path = join(self.tmp_dir, name)
		opener = gzip.open if compress else open

		with opener(path, 'wt') as fh:
			fh.write(text)

		return path

	def test_fasta_with_leading_comments(self):
		for compress in (False, True):
			fa = pyfastx.Fasta(self.write('comments{}.fa'.format(compress), COMMENTED, compress))
			self.assertEqual(list(fa.keys()), ['chr1', 'chr2'])
			self.assertEqual(fa['chr1'].seq, 'ACGTACGTACGT')
			self.assertEqual(fa.fetch('chr1', (3, 10)), 'GTACGTAC')
			self.assertEqual(fa['chr2'].seq, 'TTTTGGGG')
			del fa

	def test_fastx_with_leading_comments(self):
		path = self.write('comments.fa', COMMENTED)
		self.assertEqual(list(pyfastx.Fastx(path)), [('chr1', 'ACGTACGTACGT'), ('chr2', 'TTTTGGGG')])

	def test_comments_before_fastq_are_rejected(self):
		path = self.write('comments.fq', ";comment\n@read1\nACGT\n+\nIIII\n")

		with self.assertRaises(RuntimeError):
			pyfastx.Fastx(path)

	def test_other_text_before_the_first_record_is_rejected(self):
		path = self.write('text.fa', "not a comment\n>chr1\nACGT\n")

		with self.assertRaises(RuntimeError):
			pyfastx.Fasta(path)


if __name__ == '__main__':
	unittest.main()
