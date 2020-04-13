import os
import re
import sys
import gzip
import math
import random
import pyfastx
import argparse

def fastx_format_check(infile):
	if pyfastx.gzip_check(infile):
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
	long_cols = [0] * len(table[0])

	for row in table:
		for idx, col in enumerate(row):
			l = len(str(col))
			if l > long_cols[idx]:
				long_cols[idx] = l

	for row in table:
		row = ['{:<{}}'.format(col, long_cols[idx]) if idx==0 else
		'{:>{}}'.format(col, long_cols[idx]) for idx, col in enumerate(row)]

		print("\t".join(row))

def fastx_build(args):
	for infile in args.fastx:
		if infile.endswith('.fxi'):
			continue

		fastx_type = fastx_format_check(infile)
		
		if fastx_type == 'fasta':
			_ = pyfastx.Fasta(infile, full_index=args.full)

		elif fastx_type == 'fastq':
			_ = pyfastx.Fastq(infile, full_index=args.full)

def fastx_info(args):
	farows = [["fileName", "seqType", "seqCounts", "totalBases", "GC%",
				"meanLen", "medianLen", "maxLen", "minLen", "N50", "L50"]]
	fqrows = [["fileName", "readCounts", "totalBases", "GC%", "qualityEncodingSystem"]]
	
	for infile in args.fastx:
		if infile.endswith('.fxi'):
			continue

		fastx_type = fastx_format_check(infile)

		if fastx_type == 'fasta':
			fa = pyfastx.Fasta(infile, full_index=True)
			row = [os.path.basename(infile), fa.type, len(fa), fa.size, round(fa.gc_content, 3),
				round(fa.mean,2), round(fa.median,2), len(fa.longest), len(fa.shortest)]
			row.extend(fa.nl())
			farows.append(row)

		elif fastx_type == 'fastq':
			fq = pyfastx.Fastq(infile, full_index=True)
			row = [os.path.basename(infile), len(fq), fq.size, round(fq.gc_content,3), ",".join(fq.encoding_type)]
			fqrows.append(row)

	if len(farows) > 1:
		print_table(farows)

	if len(fqrows) > 1:
		if len(farows) > 1:
			print()
		print_table(fqrows)

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

def fasta_extract(args, ids):
	fa = pyfastx.Fasta(args.fastx)

	if args.outfile:
		fw = open(args.outfile, 'w')
	else:
		fw = sys.stdout

	for _id in ids:
		fw.write(fa[_id].raw)

	if args.outfile:
		fw.close()
	else:
		fw.flush()

def fastq_extract(args, ids):
	fq = pyfastx.Fastq(args.fastx)
	
	if args.outfile:
		fw = open(args.outfile, 'w')
	else:
		fw = sys.stdout

	if args.outfa:
		for _id in ids:
			r = fq[_id]
			fw.write("{}\n{}\n".format(r.name, r.seq))

	else:
		for _id in ids:
			fw.write(fq[_id].raw)

	if args.outfile:
		fw.close()
	else:
		fw.flush()

def fastx_extract(args):
	fastx_type = fastx_format_check(args.fastx)

	ids = None

	if args.ids:
		if re.match(r'\d+\-\d+', args.ids):
			start, end = args.ids.split('-')
			ids = range(int(start)-1, int(end))

		elif re.match(r'\d+,\d+', args.ids):
			ids = map(lambda x: int(x)-1, args.ids.split(','))

		else:
			ids = [int(args.ids)-1]

	elif args.names:
		if os.path.isfile(args.names):
			with open(args.names) as fh:
				ids = [line.strip() for line in fh]

		elif ',' in args.names:
			ids = args.names.split(',')

		else:
			ids = [args.names]

	if not ids:
		raise Exception("no ids or names input")

	if fastx_type == 'fasta':
		fasta_extract(args, ids)

	elif fastx_type == 'fastq':
		fastq_extract(args, ids)

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

	#build index command
	parser_build = subparsers.add_parser('build',
		help = "build index for FASTA or FASTQ file"
	)
	parser_build.set_defaults(func=fastx_build)
	parser_build.add_argument('-f', '--full',
		help = "build full index, base composition will be calculated",
		action = 'store_true'
	)
	parser_build.add_argument('fastx',
		help = "fasta or fastq file, gzip support",
		nargs = '+'
	)

	#statistics command
	parser_info = subparsers.add_parser('info',
		help = "show detailed statistics information of FASTA/Q file"
	)
	parser_info.set_defaults(func=fastx_info)
	parser_info.add_argument('fastx',
		help = "fasta or fastq file, gzip support",
		nargs = '+'
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
		help = "convert fastq file to fasta file"
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
		help = "get subseqence from fasta file by id or name with region"
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

	#extract sequences
	parser_extract = subparsers.add_parser('extract',
		help = "extract sequences or reads from fasta or fastq file"
	)
	parser_extract.set_defaults(func=fastx_extract)
	extract_group = parser_extract.add_mutually_exclusive_group(
		required = True
	)
	extract_group.add_argument('--ids',
		metavar = 'int or str',
		help = ("extract sequences by id number, the value can be one integer to get one sequence, "
				"a range (e.g. 5-10) or a comma seperated list (e.g. 3,5,8) to get multiple sequences")
	)
	extract_group.add_argument('--names',
		metavar = 'str',
		help = ("extract sequences by name, the value can be one name to get one sequence, "
				"a comma seperated list (e.g. seq1,seq5,seq9) or a file contains names "
				"(one name per line) to get multiple sequences")
	)
	parser_extract.add_argument('--outfas',
		help = "output fasta format when input file is fastq format, default output fastq format",
		action = 'store_true'
	)
	parser_extract.add_argument('-o', '--outfile',
		metavar = 'str',
		help = "output file, default: output to stdout"
	)
	parser_extract.add_argument('fastx',
		help = "fasta or fastq file, gzip support"
	)

	args = parser.parse_args()

	if hasattr(args, 'func'):
		args.func(args)
	else:
		parser.print_help()

if __name__ == '__main__':
	main()
	