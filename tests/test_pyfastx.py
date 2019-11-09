#!/usr/bin/env python
# coding: utf-8
from __future__ import division

import os
import gzip
import random
import pyfastx
import pyfaidx
import unittest
import statistics

#os.chdir(os.path.dirname(__file__))

gzip_fasta = 'tests/data/test.fa.gz'
flat_fasta = 'tests/data/test.fa'
gzip_fastq = 'tests/data/test.fq.gz'

class FastaTest(unittest.TestCase):
	def setUp(self):
		self.fastx = pyfastx.Fasta(gzip_fasta)

		#reload index
		self.fastx = pyfastx.Fasta(gzip_fasta)

		self.faidx = pyfaidx.Fasta(flat_fasta, sequence_always_upper=True)
		self.fastq = pyfastx.Fastq(gzip_fastq)
		self.count = len(self.fastx)

		self.reads = {}
		self.bases = {'A': 0, 'T': 0, 'G': 0, 'C':0, 'N':0}
		i = 0
		c = -1
		with gzip.open(gzip_fastq, 'rt') as fh:
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
		if os.path.exists('tests/data/test.fa.gz.fxi'):
			os.remove('tests/data/test.fa.gz.fxi')

		if os.path.exists('tests/data/test.fa.fai'):
			os.remove('tests/data/test.fa.fai')

		if os.path.exists('tests/data/test.fq.gz.fxi'):
			os.remove('tests/data/test.fq.gz.fxi')

	def get_random_index(self):
		return random.randint(0, self.count-1)

	def get_random_read(self):
		return random.randint(0, len(self.fastq)-1)

	def test_module(self):
		# gzip check test
		self.assertEqual(pyfastx.gzip_check(gzip_fasta), self.fastx.is_gzip)

		# version test
		with open('src/version.h') as fh:
			version = fh.read().split()[2].strip('"')
			self.assertEqual(version, pyfastx.version())

	def test_fasta(self):
		#fasta format
		self.assertEqual(self.fastx.type, 'DNA')

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

		self.assertEqual(longest, self.fastx.longest)
		self.assertEqual(shortest, self.fastx.shortest)

		#test contains
		idx = self.get_random_index()
		name = self.faidx[0].name
		self.assertTrue(name in self.fastx)

	def test_iter(self):
		for name, result in self.fastx:
			expect = self.faidx[name][:].seq
			self.assertEqual(expect, result)

	def test_stat(self):
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
		expect = round(statistics.mean(lens), 3)
		result = round(self.fastx.mean, 3)
		self.assertEqual(expect, result)

		#test median length
		expect = statistics.median(lens)
		result = self.fastx.median
		self.assertEqual(expect, result)

		#test count squence
		expect = 0
		for l in lens:
			if l >= 200:
				expect += 1
		result = self.fastx.count(200)
		self.assertEqual(expect, result)

	def test_keys(self):
		expect = list(self.faidx.keys())
		result = list(self.fastx.keys())

		self.assertEqual(expect, result)

		#id counts
		ids = self.fastx.keys()
		self.assertEqual(len(ids), len(expect))

		#get id from identifier class
		self.assertEqual(ids[0], expect[0])
		self.assertEqual(ids[-1], expect[-1])

		#check contains
		idx = self.get_random_index()
		name = self.faidx[idx].name

		self.assertTrue(name in ids)

		#sor by id
		expect = [seq.name for seq in self.faidx]
		expect.reverse()
		result = [name for name in self.fastx.keys().sort('id', reverse=True)]
		self.assertEqual(expect, result)

		#sort by name
		expect = sorted([seq.name for seq in self.faidx])
		result = [name for name in self.fastx.keys().sort('name')]
		self.assertEqual(expect, result)

		#sort by length
		lens = [(seq.name, len(seq)) for seq in self.faidx]
		expect = [it[0] for it in sorted(lens, key=lambda x: x[1])]
		result = [name for name in self.fastx.keys().sort('length')]
		self.assertEqual(expect, result)

	def test_seq_by_index(self):
		#test get seq by index
		idx = self.get_random_index()
		expect = self.faidx[idx][:]
		result = self.fastx[idx]

		self.assertEqual(expect.name, result.name)
		self.assertEqual(expect.seq, result.seq)

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

	def test_seq_slice(self):
		idx = self.get_random_index()
		expect = self.faidx[idx]
		result = self.fastx[idx]

		self.assertEqual(expect[5:10].seq, result[5:10].seq)
		#self.assertEqual(expect[20:].seq, result[20:].seq)
		expect = expect[20:].seq
		result = result[20:].seq
		if len(expect) > len(result):
			expect = expect[0:-1]
		
		self.assertEqual(expect, result)

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

	def test_seq_iter(self):
		idx = self.get_random_index()
		fai_seq = self.faidx[idx]
		fxi_seq = self.fastx[idx]

		# test read seq line by line
		for expect, result in zip(fai_seq, fxi_seq):
			self.assertEqual(str(expect), result)

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

	def test_get_seq(self):
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

		# test length
		self.assertEqual(len(result), len(expect[1]))

		# test name
		self.assertEqual(result.name, expect[0])

		# test str
		self.assertEqual(str(result), expect[1])

		# test seq
		self.assertEqual(result.seq, expect[1])

		# test quality
		self.assertEqual(result.qual, expect[2])

		# test quality integer
		self.assertEqual(result.quali, [ord(b)-33 for b in expect[2]])

		result = self.fastq[expect[0]]

		# test subscript
		self.assertEqual(result.seq, expect[1])

		# test contain
		self.assertTrue(result.name in self.fastq)

		# test read iter
		i = -1
		for read in self.fastq:
			i += 1
			self.assertEqual(read.name, self.reads[i][0])
			self.assertEqual(read.seq, self.reads[i][1])
			self.assertEqual(read.qual, self.reads[i][2])
