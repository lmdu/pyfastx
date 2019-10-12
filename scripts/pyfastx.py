import os
import sys
import gzip
import math
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

def fasta_split(args):
	fa = pyfastx.Fasta(args.fastx)

	if args.seq_count:
		parts_num = math.ceil(len(fa)/args.seq_count)
	else:
		parts_num = args.file_num

	name, suffix = os.path.splitext(os.path.basename(args.fastx))
	digit = len(str(parts_num))
	lens = [0] * parts_num
	
	if args.seq_count:
		seqs = [0] * parts_num

	fhs = []

	for i in range(1, parts_num+1):
		subfile = "{}.{}{}".format(name, str(i).zfill(digit), suffix)

		if args.out_dir != '.':
			subfile = os.path.join(args.out_dir, subfile)

		if fa.gzip:
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

	name, suffix = os.path.splitext(os.path.basename(args.fastx))
	digit = len(str(parts_num))

	seq_write = 0
	fh = None
	file_num = 0

	for read in fq:
		if seq_write == 0:
			file_num += 1
			subfile = "{}.{}{}".format(name, str(file_num).zfill(digit), suffix)

			if args.out_dir != '.':
				subfile = os.path.join(args.out_dir, subfile)

			if fq.gzip:
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


if __name__ == '__main__':
	parser = argparse.ArgumentParser(
		prog = 'pyfastx',
		usage = "pyfastx COMMAND [OPTIONS]",
		description = "A tool for FASTA/Q file manipulation",
		epilog = "Contact:\nLianming Du (dulianming@cdu.edu.cn)",
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

	parser_split = subparsers.add_parser('split', 
		help = "Split fasta file into multiple files"
	)
	parser_split.set_defaults(func=fastx_split)

	group = parser_split.add_mutually_exclusive_group(required=True)

	group.add_argument('-n', '--file_num',
		dest = 'file_num',
		type = int,
		default = 0,
		help = "split a fasta or fastq file into N new files with even size"
	)

	group.add_argument('-c', '--seq_count',
		dest = 'seq_count',
		type = int,
		default = 0,
		help = "split a fasta or fastq file into multiple files with the same sequence counts"
	)

	parser_split.add_argument('-o', '--out_dir',
		dest = 'out_dir',
		help = "output directory, default is current folder",
		default = '.'
	)

	parser_split.add_argument('-g', '--gzip_compress',
		dest = 'out_gzip',
		action = 'store_true',
		default = False,
		help = 'use gzip to compress output files',
	)

	parser_split.add_argument('fastx',
		help = 'input fasta or fastq file, gzip compressed support'
	)


	args = parser.parse_args()

	try:
		args.func(args)

	except AttributeError:
		parser.print_help()


