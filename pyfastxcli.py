import os
import re
import sys
import gzip
import math
import random
import pyfastx
import argparse

def fastx_format_check(infile):
	if (pyfastx.gzip_check(infile)):
		fp = gzip.open(infile, 'rt')

	else:
		fp = open(infile)
	
	for line in fp:
		if line.strip():
			break
	
	fp.close()

	if line[0] == '>':
		return 'fasta'

	elif line[0] == '@':
		return 'fastq'

	else:
		raise Exception("Input file %s is not fasta or fastq", infile)

#find max and min value index
def max_index(lst):
	return max(range(len(lst)), key=lst.__getitem__)

def min_index(lst):
	return min(range(len(lst)), key=lst.__getitem__)

def is_glob(pattern):
	if re.search(r'[\*\?]', pattern) or re.search(r'\[\!?.*\]', pattern):
		return True
	else:
		return False

def print_table(table):
	width = [0] * len(table[0])

	for row in table:
		for i, col in enumerate(row):
			l = len(str(col))

			if l > width[i]:
				width[i] = l

	format_row = "\t".join(['{:%s}' % j for j in width])

	for row in table:
		print(format_row.format(row))

def fastx_info(args):
	fastx_type = fastx_format_check(args.fastx)

	if fastx_type == 'fasta':
		fa = pyfastx.Fasta(args.fastx)
		comp = fa.composition
		print("Sequence counts: {}".format(len(fa)))
		print("Total bases: {}".format(fa.size))
		print("GC content: {:.2f}%".format(fa.gc_content))
		for b in comp:
			print("{} counts: {}".format(b, comp[b]))
		print("Mean length: {:.2f}".format(fa.mean))
		print("Median length: {:.2f}".format(fa.median))
		print("Max length: {}".format(len(fa.longest)))
		print("Min length: {}".format(len(fa.shortest)))
		print("N50, L50: {}, {}".format(*fa.nl()))
		print("length >= 1000: {}".format(fa.count(1000)))

	elif fastx_type == 'fastq':
		fq = pyfastx.Fastq(args.fastx)
		comp = fq.composition
		print("Read counts: {}".format(len(fq)))
		print("Total bases: {}".format(fq.size))
		print("GC content: {:.2f}%".format(fq.gc_content))
		for b in comp:
			print("{} counts: {}".format(b, comp[b]))
		print("Quality encoding system maybe: {}".format(", ".join(fq.encoding_type)))

def fasta_split(args):
	fa = pyfastx.Fasta(args.fastx)

	if args.seq_count:
		parts_num = math.ceil(len(fa)/args.seq_count)
	else:
		parts_num = args.file_num

	name, suffix1 = os.path.splitext(os.path.basename(args.fastx))

	if fa.is_gzip:
		name, suffix2 = os.path.splitext(name)

	digit = len(str(parts_num))
	lens = [0] * parts_num
	
	if args.seq_count:
		seqs = [0] * parts_num

	fhs = []

	for i in range(1, parts_num+1):
		if fa.is_gzip:
			subfile = "{}.{}{}{}".format(name, str(i).zfill(digit), suffix2, suffix1)
		else:
			subfile = "{}.{}{}".format(name, str(i).zfill(digit), suffix1)

		if args.outdir is not None:
			subfile = os.path.join(args.outdir, subfile)

		if fa.is_gzip:
			fh = gzip.open(subfile, 'wt')
		else:
			fh = open(subfile, 'w')

		fhs.append(fh)

	ids = fa.keys()
	for chrom in ids.sort('length', reverse=True):
		idx = min_index(lens)
		fhs[idx].write(">%s\n%s\n" % (chrom, fa[chrom].seq))
		lens[idx] += len(fa[chrom])

		if args.seq_count:
			seqs[idx] += 1

			if seqs[idx] == args.seq_count:
				lens[idx] = fa.size

	for fh in fhs:
		fh.close()

def fastq_split(args):
	fq = pyfastx.Fastq(args.fastx)

	if args.file_num:
		seqs_num = math.ceil(len(fq)/args.file_num)
		parts_num = args.file_num
	else:
		seqs_num = args.seq_count
		parts_num = math.ceil(len(fq)/seqs_num)

	name, suffix1 = os.path.splitext(os.path.basename(args.fastx))

	if fq.is_gzip:
		name, suffix2 = os.path.splitext(name)

	digit = len(str(parts_num))

	seq_write = 0
	fh = None
	file_num = 0

	for read in fq:
		if seq_write == 0:
			file_num += 1

			if fq.is_gzip:
				subfile = "{}.{}{}{}".format(name, str(file_num).zfill(digit), suffix2, suffix1)
			else:
				subfile = "{}.{}{}".format(name, str(file_num).zfill(digit), suffix1)

			if args.outdir is not None:
				subfile = os.path.join(args.outdir, subfile)

			if fq.is_gzip:
				fh = gzip.open(subfile, 'wt')
			else:
				fh = open(subfile, 'w')

		fh.write("@{}\n{}\n+\n{}\n".format(read.name, read.seq, read.qual))
		seq_write += 1

		if seq_write == seqs_num:
			fh.close()
			seq_write = 0
	
	fh.close()

def fastx_split(args):
	fastx_type = fastx_format_check(args.fastx)

	if fastx_type == 'fasta':
		fasta_split(args)

	elif fastx_type == 'fastq':
		fastq_split(args)

def fastx_fq2fa(args):
	fq = pyfastx.Fastq(args.fastx)

	if args.outfile:
		fh = open(args.outfile, 'w')
	else:
		fh = sys.stdout

	for read in fq:
		fh.write(">{}\n{}\n".format(read.name, read.seq))

	if args.outfile:
		fh.close()
	else:
		fh.flush()

def fastx_subseq(args):
	fa = pyfastx.Fasta(args.fastx)

	if args.chr is not None:
		if args.chr not in fa:
			raise RuntimeError("no sequence named {} in fasta file".format(args.chr))

		subseq = fa[args.chr]

	else:
		if args.id <= 0:
			raise RuntimeError("sequence id must be a integer between 1 and {}".format(len(fa)))

		subseq = fa[args.id]

	if args.region:
		start, end = args.region.split(':')
		if start:
			start = int(start) - 1
		else:
			start = 0

		if end:
			end = int(end)
		else:
			end = len(s)

		sys.stdout.write("{}\n".format(subseq[start:end].seq))
	else:
		sys.stdout.write("{}\n".format(subseq.seq))

	sys.stdout.flush()

def fasta_sample(args):
	fa = pyfastx.Fasta(args.fastx)

	if args.num is not None and args.num > 0:
		seq_num = args.num
		if seq_num > len(fa):
			seq_num = len(fa)

	elif args.prop is not None and 0 < args.prop <= 1:
		seq_num = round(len(fa)*args.prop)
		if seq_num == 0:
			raise RuntimeError("the proportion is too small")

	else:
		raise RuntimeError("specify a right number for seq number or proportion")

	selected = random.sample(range(len(fa)), k=seq_num)

	if args.outfile is None:
		fw = sys.stdout
	else:
		fw = open(args.outfile, 'w')

	for idx in selected:
		s = fa[idx]
		fw.write(">{}\n{}\n".format(s.name, s.seq))

	if args.outfile is None:
		fw.flush()
	else:
		fw.close()

def fastq_sample(args):
	fq = pyfastx.Fastq(args.fastx)

	if args.num is not None and args.num > 0:
		seq_num = args.num
		if seq_num > len(fq):
			seq_num = len(fq)

	elif args.num is not None and 0 < args.prop <= 1:
		seq_num = round(len(fq)*args.prop)
		if seq_num == 0:
			raise RuntimeError("the proportion is too small")

	else:
		raise RuntimeError("specify a right number for seq number or proportion")

	selected = random.sample(range(len(fq)), k=seq_num)

	if args.outfile is None:
		fw = sys.stdout
	else:
		fw = open(args.outfile, 'w')

	for idx in selected:
		r = fq[idx]
		fw.write("@{}\n{}\n+\n{}\n".format(r.name, r.seq, r.qual))

	if args.outfile is None:
		fw.flush()
	else:
		fw.close()

def fastx_sample(args):
	fastx_type = fastx_format_check(args.fastx)

	if fastx_type == 'fasta':
		fasta_sample(args)

	elif fastx_type == 'fastq':
		fastq_sample(args)

def main():
	parser = argparse.ArgumentParser(
		prog = 'pyfastx',
		usage = "pyfastx COMMAND [OPTIONS]",
		description = "A command line tool for FASTA/Q file manipulation",
		formatter_class = argparse.RawDescriptionHelpFormatter
	)

	parser.add_argument('-v', '--version',
		action = 'version',
		version = "%(prog)s version {}".format(pyfastx.version())
	)

	subparsers = parser.add_subparsers(
		title = 'Commands',
		prog = 'pyfastx',
		metavar = ''
	)

	#statistics command
	parser_info = subparsers.add_parser('info',
		help = "show detailed statistics information of FASTA/Q file"
	)
	parser_info.set_defaults(func=fastx_info)
	parser_info.add_argument('fastx',
		help = "fasta or fastq file, gzip support"
	)

	#split command
	parser_split = subparsers.add_parser('split',
		help = "split fasta file into multiple files"
	)
	parser_split.set_defaults(func=fastx_split)
	split_group = parser_split.add_mutually_exclusive_group(
		required=True
	)
	split_group.add_argument('-n',
		dest = 'file_num',
		type = int,
		metavar = 'int',
		help = "split a fa/q file into N new files with even size"
	)
	split_group.add_argument('-c',
		dest = 'seq_count',
		type = int,
		metavar = 'int',
		help = "split a fa/q file into multiple files with the same sequence counts"
	)
	parser_split.add_argument('-o', '--outdir',
		dest = 'outdir',
		help = "output directory, default is current folder",
		metavar = 'str'
	)
	parser_split.add_argument('fastx',
		help = 'fasta or fastq file, gzip support'
	)

	#convert fastq to fasta command
	parser_fq2fa = subparsers.add_parser('fq2fa',
		help = "Convert fastq file to fasta file"
	)
	parser_fq2fa.set_defaults(func=fastx_fq2fa)
	parser_fq2fa.add_argument('-o', '--outfile',
		dest = 'outfile',
		metavar = 'str',
		help = "output file, default: output to stdout"
	)

	parser_fq2fa.add_argument('fastx',
		help = "input fastq file, gzip support"
	)

	#get subseq from fasta
	parser_subseq = subparsers.add_parser('subseq',
		help = "Get subseqence from fasta file by id or name with region"
	)
	parser_subseq.set_defaults(func=fastx_subseq)

	subseq_group = parser_subseq.add_mutually_exclusive_group(
		required = True
	)

	subseq_group.add_argument('--id',
		help = "sequence id number in fasta file",
		type = int,
		metavar = 'int'
	)
	subseq_group.add_argument('--chr',
		help = "sequence name",
		metavar = 'str'
	)
	parser_subseq.add_argument('-r', '--region',
		help = "one-based slice region, e.g. 10:20",
		metavar = 'str'
	)
	parser_subseq.add_argument('fastx',
		help = "input fasta file, gzip support"
	)

	parser_sample = subparsers.add_parser('sample',
		help = "randomly sample sequences from fasta or fastq file"
	)
	parser_sample.set_defaults(func=fastx_sample)
	sample_group = parser_sample.add_mutually_exclusive_group(
		required = True
	)
	sample_group.add_argument('-n',
		dest = 'num',
		help = "number of sequences to be sampled",
		type = int,
		metavar = 'int'
	)
	sample_group.add_argument('-p',
		dest = 'prop',
		help = "proportion of sequences to be sampled, 0~1",
		type = float,
		metavar = 'float'
	)
	parser_sample.add_argument('-o', '--outfile',
		metavar = 'str',
		help = "output file, default: output to stdout"
	)
	parser_sample.add_argument('fastx',
		help = "fasta or fastq file, gzip support"
	)

	args = parser.parse_args()

	if hasattr(args, 'func'):
		args.func(args)
	else:
		parser.print_help()

if __name__ == '__main__':
	main()
	