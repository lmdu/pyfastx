import os
import random
import pyfastx
import pyfaidx
import unittest

join = os.path.join
root_dir = join(os.path.dirname(os.path.abspath(__file__)), '..')
data_dir = join(os.path.dirname(os.path.abspath(__file__)), 'data')

gzip_fasta = join(data_dir, 'test.fa.gz')
flat_fasta = join(data_dir, 'test.fa')
rna_fasta = join(data_dir, 'rna.fa')
protein_fasta = join(data_dir, 'protein.fa')


class FastaTest(unittest.TestCase):
	def setUp(self):
		self.fastx = pyfastx.Fasta(gzip_fasta)

		self.fasta = pyfastx.Fasta(flat_fasta)

		self.faidx = pyfaidx.Fasta(flat_fasta, sequence_always_upper=True)

		self.count = len(self.fastx)

	def tearDown(self):
		#fix permission error on Windows
		del self.fastx
		del self.fasta
		del self.faidx

		if os.path.exists('{}.fxi'.format(gzip_fasta)):
			os.remove('{}.fxi'.format(gzip_fasta))

		if os.path.exists('{}.fxi'.format(flat_fasta)):
			os.remove('{}.fxi'.format(flat_fasta))

		if os.path.exists('{}.fai'.format(flat_fasta)):
			os.remove('{}.fai'.format(flat_fasta))

		if os.path.exists('{}.fxi'.format(rna_fasta)):
			os.remove('{}.fxi'.format(rna_fasta))

		if os.path.exists('{}.fxi'.format(protein_fasta)):
			os.remove('{}.fxi'.format(protein_fasta))

	def get_random_index(self):
		return random.randint(0, self.count-1)

	def test_module(self):
		# gzip check test
		self.assertEqual(pyfastx.gzip_check(gzip_fasta), self.fastx.is_gzip)

		# version test
		with open(join(root_dir, 'src', 'version.h')) as fh:
			version = fh.read().split()[2].strip('"')
			self.assertEqual(version, pyfastx.version())

		print(pyfastx.version(debug=True))

	def test_build(self):
		self.fastx = pyfastx.Fasta(gzip_fasta, build_index=False)

		if os.path.exists('{}.fxi'.format(gzip_fasta)):
			os.remove('{}.fxi'.format(gzip_fasta))

		self.fastx.build_index()

	def test_fasta(self):
		#test gzip
		self.assertFalse(self.fasta.is_gzip)

		#seq counts
		self.assertEqual(len(self.fastx), len(self.faidx.keys()))

		#seq length
		expect_size = sum(len(s) for s in self.faidx)
		self.assertEqual(self.fastx.size, expect_size)

		#test composition
		expect = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
		for s in self.faidx:
			expect['A'] += s[:].seq.count('A')
			expect['C'] += s[:].seq.count('C')
			expect['G'] += s[:].seq.count('G')
			expect['T'] += s[:].seq.count('T')
		self.assertEqual(self.fastx.composition, expect)

		#test GC content
		expect_gc = (expect['G']+expect['C'])/sum(expect.values())*100
		self.assertEqual(round(self.fastx.gc_content, 3), round(expect_gc, 3))

		#test GC skew
		expect_skew = (expect['G']-expect['C'])/(expect['G']+expect['C'])
		self.assertEqual(round(self.fastx.gc_skew, 3), round(expect_skew, 3))

		#test longest and shortest sequence
		longest = (None, 0)
		shortest = (None, expect_size)
		for seq in self.faidx:
			l = len(seq)
			if l > longest[1]:
				longest = (seq.name, l)

			if l < shortest[1]:
				shortest = (seq.name, l)

		long_seq = self.fastx.longest
		short_seq = self.fastx.shortest

		self.assertEqual(longest, (long_seq.name, len(long_seq)))
		self.assertEqual(shortest, (short_seq.name, len(short_seq)))

		#test contains
		idx = self.get_random_index()
		name = self.faidx[idx].name
		self.assertTrue(name in self.fastx)

	#test repr
	def test_repr(self):
		expect = "<Fasta> {} contains {} sequences".format(gzip_fasta, self.count)
		result = repr(self.fastx)
		self.assertEqual(expect, result)

		#without build index
		fa = pyfastx.Fasta(flat_fasta, build_index=False)
		expect = "<Fasta> {}".format(flat_fasta)
		result = repr(fa)
		self.assertEqual(expect, result)

	def test_seq_type(self):
		#test dna format
		self.assertEqual(self.fastx.type, 'DNA')

		#test rna format
		rna = pyfastx.Fasta(rna_fasta)
		self.assertEqual(rna.type, "RNA")

		#test protein format
		prot = pyfastx.Fasta(protein_fasta)
		self.assertEqual(prot.type, "protein")

	def test_iter_object(self):
		for seq in self.fastx:
			expect = self.faidx[seq.name][:].seq
			self.assertEqual(expect, seq.seq)

		for seq in pyfastx.Fasta(flat_fasta, uppercase=True):
			expect = self.faidx[seq.name][:].seq
			self.assertEqual(expect, seq.seq)

		#test reference of sequence made from loop
		for seq in self.fastx:
			break

		expect = self.faidx[seq.name][:].seq
		self.assertEqual(expect, seq.seq)
		self.assertEqual(expect, seq.seq)

	def test_iter_tuple(self):
		fa = pyfastx.Fasta(gzip_fasta, build_index=False)
		
		for name, seq in fa:
			expect = str(self.faidx[name])
			self.assertEqual(expect, seq)

	def test_iter_upper(self):
		fa = pyfastx.Fasta(flat_fasta, build_index=False, uppercase=True)

		for name, seq in fa:
			self.assertEqual(seq, self.faidx[name.split()[0]][:].seq)

	def test_iter_full_name(self):
		fa = pyfastx.Fasta(flat_fasta, build_index=False, full_name=True)

		for name, _ in fa:
			self.assertTrue(name, self.fastx[name.split()[0]].description)

	def test_iter_upper_full_name(self):
		fa = pyfastx.Fasta(flat_fasta, build_index=False, uppercase=True, full_name=True)

		for name, seq in fa:
			self.assertEqual(name, self.fastx[name.split()[0]].description)
			self.assertEqual(seq, self.fastx[name.split()[0]].seq)

	def test_key_func(self):
		del self.fastx

		#remove previously created index file
		if os.path.exists("{}.fxi".format(gzip_fasta)):
			os.remove("{}.fxi".format(gzip_fasta))

		self.fastx = pyfastx.Fasta(gzip_fasta, key_func=lambda x: x.split()[1])
		idx = self.get_random_index()
		self.assertEqual(self.fastx[idx].name, self.fastx[idx].description.split()[1])

	def test_statistics(self):
		lens = sorted([len(seq) for seq in self.faidx], reverse=True)
		half = sum(lens)/2
		tlen = 0
		l50 = 0
		for n50 in lens:
			l50 += 1
			tlen += n50

			if tlen >= half:
				break

		self.assertEqual(self.fastx.nl(50), (n50, l50))

		#test mean length
		expect = round(sum(lens)/len(lens), 3)
		result = round(self.fastx.mean, 3)
		self.assertEqual(expect, result)

		#test median length
		lens = sorted(lens)
		expect = lens[105]

		result = self.fastx.median
		self.assertEqual(expect, result)

		#test count squence
		expect = 0
		for l in lens:
			if l >= 200:
				expect += 1
		result = self.fastx.count(200)
		self.assertEqual(expect, result)

	def test_seq_fetch(self):
		idx = self.get_random_index()
		name = list(self.faidx.keys())[idx]
		l = len(self.fastx[idx])

		#test one interval
		a = int(l/2)
		interval = (random.randint(1, a), random.randint(a+1, l))

		expect = str(self.faidx[name])[interval[0]-1:interval[1]]
		result = self.fastx.fetch(name, interval)

		self.assertEqual(expect, result)

		#test multiple intervals
		intervals = []
		intervals.append((random.randint(1, int(a/2)), random.randint(int(a/2)+1, a)))
		intervals.append((random.randint(a+1, int((a+l)/2)), random.randint(int((a+l)/2)+1, l)))

		expect = "".join([str(self.faidx[name])[s-1:e] for s, e in intervals])
		result = self.fastx.fetch(name, intervals)

		self.assertEqual(expect, result)

	def test_seq_flank(self):
		idx = self.get_random_index()
		name = list(self.faidx.keys())[idx]
		l = len(self.fastx[idx])

		a = int(l/2)
		start = random.randint(1, a)
		end = random.randint(a+1, l)
		flen = 20
		left, right = self.fastx.flank(name, start, end, flen)
		s = start - flen - 1
		if s < 0:
			s = 0

		self.assertEqual(str(self.faidx[name])[s:start-1], left)
		self.assertEqual(str(self.faidx[name])[end:end+flen], right)

		left, right = self.fastx.flank(name, start, end, flank_length=flen, use_cache=True)
		self.assertEqual(str(self.faidx[name])[s:start-1], left)
		self.assertEqual(str(self.faidx[name])[end:end+flen], right)

		left, right = self.fastx.flank(name, 1, len(self.faidx[name]))
		self.assertEqual('', left)
		self.assertEqual('', right)

	def test_no_upper(self):
		fa = pyfastx.Fasta(flat_fasta, uppercase=False)
		self.assertEqual(self.fastx[0].seq, fa[0].seq)

	def test_exception(self):
		with self.assertRaises(TypeError):
			pyfastx.Fasta(flat_fasta, key_func=1)

		with self.assertRaises(FileExistsError):
			pyfastx.Fasta('a_file_not_exists')

		with self.assertRaises(ValueError):
			self.fastx.fetch('seq1', {'a':1})

		with self.assertRaises(NameError):
			self.fastx.fetch('seq1', (1,10))

		with self.assertRaises(ValueError):
			self.fastx.fetch(self.fastx[0].name, (1,10,20))

		with self.assertRaises(ValueError):
			self.fastx.fetch(self.fastx[0].name, (20, 10))

		with self.assertRaises(ValueError):
			self.fastx.fetch(self.fastx[0].name, [20, 10])

		with self.assertRaises(IndexError):
			_ = self.fastx[self.count]

		with self.assertRaises(KeyError):
			_ = self.fastx[list()]

		with self.assertRaises(ValueError):
			self.fastx.nl(101)

		with self.assertRaises(RuntimeError):
			non_fa = 'non.fa'
			with open(non_fa, 'w') as fw:
				fw.write('abc')

			pyfastx.Fasta(non_fa)

			os.remove(non_fa)

if __name__ == '__main__':
	unittest.main()
