import os
import random
import pyfastx
import pyfaidx
import unittest

join = os.path.join
data_dir = join(os.path.dirname(os.path.abspath(__file__)), 'data')

gzip_fasta = join(data_dir, 'test.fa.gz')
flat_fasta = join(data_dir, 'test.fa')

class SequenceTest(unittest.TestCase):
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

	def test_seq_by_index(self):
		#test get seq by index
		idx = self.get_random_index()
		expect = self.faidx[idx][:]
		result = self.fastx[idx]

		self.assertEqual(expect.name, result.name)
		self.assertEqual(expect.seq, result.seq)

		#test subseq
		self.assertEqual(expect[0:10].seq, result[0:10].seq)

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
		self.assertEqual(expect.seq, result.seq)

	def test_seq_reverse_complement(self):
		idx = self.get_random_index()
		expect = self.faidx[idx][:]
		result = self.fastx[idx]

		self.assertEqual(str(expect.reverse), result.reverse)
		self.assertEqual(str(expect.complement), result.complement)
		self.assertEqual(str(-expect), result.antisense)

	def test_seq_raw(self):
		idx = self.get_random_index()

		seq = self.fastx[idx]

		lines = []
		with open(flat_fasta) as fh:
			for line in fh:
				if line.startswith('>{}'.format(seq.name)):
					lines.append(line)
					continue

				if lines:
					if line[0] == '>':
						break

					lines.append(line)

		self.assertEqual(''.join(lines).replace('\n','\r\n'), seq.raw)

	def test_seq_slice(self):
		idx = self.get_random_index()
		expect = str(self.faidx[idx])
		result = self.fastx[idx]
		flatseq = self.fasta[idx]

		#test gzip subseq
		self.assertEqual(expect[5:10], result[5:10].seq)

		#test flat subseq
		self.assertEqual(expect[5:10], flatseq[5:10].seq)
		self.assertEqual(expect[20:], result[20:].seq)

		#test two level slice
		self.assertEqual(expect[10:100][:20], result[10:100][:20].seq)

		#test sequence index
		pos = random.randint(0, len(expect) - 1)
		self.assertEqual(expect[pos], result[pos])

		pos = random.randint(1, len(expect)) * -1
		self.assertEqual(expect[pos], result[pos])

		del flatseq

	def test_seq_content(self):
		idx = self.get_random_index()
		result = self.fastx[idx]
		expect = self.faidx[idx]

		content = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
		content['A'] += expect[:].seq.count('A')
		content['C'] += expect[:].seq.count('C')
		content['G'] += expect[:].seq.count('G')
		content['T'] += expect[:].seq.count('T')

		expect_gc = (content['G']+content['C'])/sum(content.values())*100

		self.assertEqual(result.composition, content)
		self.assertEqual(round(result.gc_content, 3), round(expect_gc, 3))

		#test gc skew
		expect_skew = (content['G']-content['C'])/(content['G']+content['C'])
		self.assertEqual(round(result.gc_skew, 3), round(expect_skew, 3))

	def test_seq_iter(self):
		idx = self.get_random_index()
		fai_seq = self.faidx[idx]
		fxi_seq = self.fastx[idx]
		fas_seq = self.fasta[idx]

		# test read seq line by line
		flatsq = [line for line in fas_seq]
		expect = [str(line) for line in fai_seq]
		result = [line for line in fxi_seq]

		self.assertEqual(expect, result)
		self.assertEqual(expect, flatsq)

		# test seq long name
		self.assertEqual(fai_seq.long_name.strip(), fxi_seq.description)

		# test seq str
		self.assertEqual(str(fai_seq), str(fxi_seq))

		# test seq contains
		s, e = sorted(random.sample(range(1, len(fai_seq)), 2))
		segment = str(fai_seq)[s-1:e]
		self.assertTrue(segment in fxi_seq)

		# test seq search
		expect = str(fai_seq).index(segment) + 1
		result = fxi_seq.search(segment)
		self.assertEqual(expect, result)

	def test_seq_repr(self):
		s = self.fastx[0]
		n = s.name
		self.assertEqual(repr(s), "<Sequence> {} with length of {}".format(s.name, len(s)))

		s = self.fastx[0][0:11]
		self.assertEqual(repr(s), "<Sequence> {} from {} to {}".format(n, s.start, s.end))

		#test name
		self.assertEqual(s.name, "{}:{}-{}".format(n, s.start, s.end))

	def test_full_compo(self):
		idx = self.get_random_index()
		fas = pyfastx.Fasta(flat_fasta, full_index=True)
		self.assertEqual(self.fastx[idx].composition, fas[idx].composition)

	def test_seq_exception(self):
		with self.assertRaises(RuntimeError):
			for line in self.fastx[0][10:20]:
				print(line)

if __name__ == '__main__':
	unittest.main()
