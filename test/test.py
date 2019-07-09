import os
import pyfastx
import sqlite3
import unittest

os.chdir(os.path.dirname(__file__))

fasta_file = 'data/test.fa.gz'
fasta_index = 'data/test.fa.gz.db'

class FastaTest(unittest.TestCase):
	def setUp(self):
		self.fasta = pyfastx.Fasta(fasta_file)
		self.conn = sqlite3.connect(fasta_index)

	def tearDown(self):
		try:
			os.remove('data/test.fa.gz.db')
		except EnvironmentError:
			pass

	def expect_val(self, sql):
		cursor = self.conn.cursor()
		cursor.execute(sql)
		row = cursor.fetchone()
		cursor.close()
		if len(row) == 1:
			return row[0]
		else:
			return row

	def test_build_index(self):
		self.assertTrue(os.path.exists(fasta_index))

	def test_fasta_attrs(self):
		fasta = self.fasta

		sql = "SELECT COUNT(*) FROM seq LIMIT 1"
		self.assertEqual(expect_val(sql), fasta.count)

		sql = "SELECT SUM(slen) FROM seq LIMIT 1"
		self.assertEqual(expect_val(sql), fasta.size)

		self.assertEqual(fasta_file, fasta.file_name)

	def test_get_seq(self):
		fasta = self.fasta

		seq = fasta[0]

		sql = "SELECT * FROM seq WHERE id=0 LIMIT 1"
		row = expect_val(sql)

		self.assertEqual(seq.name, row[1])
		self.assertEqual(len(seq.seq), row[4])



