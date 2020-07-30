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
				"avgLen", "medianLen", "maxLen", "minLen", "N50", "L50"]]
	fqrows = [["fileName", "readCounts", "totalBases", "GC%", "avgLen", "maxLen",
				"minLen", "maxQual", "minQual", "qualEncodingSystem"]]

	for infile in args.fastx:
		if infile.endswith('.fxi'):
			continue

		fastx_type = fastx_format_check(infile)

		if fastx_type == 'fasta':
			fa = pyfastx.Fasta(infile, full_index=True)
			if fa.type in ['DNA', 'RNA']:
				gc = round(fa.gc_content, 3)
			else:
				gc = '-'

			row = [os.path.basename(infile), fa.type, len(fa), fa.size, gc, round(fa.mean,3),
					round(fa.median,3), len(fa.longest), len(fa.shortest)]
			row.extend(fa.nl())
			farows.append(row)

		elif fastx_type == 'fastq':
			fq = pyfastx.Fastq(infile, full_index=True)
			row = [os.path.basename(infile), len(fq), fq.size, round(fq.gc_content,3), round(fq.avglen,3),
					fq.maxlen, fq.minlen, fq.maxqual, fq.minqual, ",".join(fq.encoding_type)]
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
	mapping = {}
	for chrom in ids.sort('length', reverse=True):
		idx = min_index(lens)
		#fhs[idx].write(fa[chrom].raw)
		mapping[chrom] = idx
		lens[idx] += len(fa[chrom])

		if args.seq_count:
			seqs[idx] += 1

			if seqs[idx] == args.seq_count:
				lens[idx] = fa.size

	for seq in fa:
		fhs[mapping[seq.name]].write(seq.raw)

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

		fh.write(read.raw)
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
	fq = pyfastx.Fastq(args.fastx, build_index=False)

	if args.out_file:
		fh = open(args.out_file, 'w')
	else:
		fh = sys.stdout

	for name, seq, _ in fq:
		fh.write(">{}\n{}\n".format(name, seq))

	if args.out_file:
		fh.close()
	else:
		fh.flush()


def fastx_subseq(args):
	fa = pyfastx.Fasta(args.fastx)

	if args.out_file:
		fw = open(args.out_file, 'w')
	else:
		fw = sys.stdout

	if args.region_file:
		with open(args.region_file) as fh:
			for line in fh:
				chrom, start, end = line.strip().split()
				start = int(start)
				end = int(end)
				seq = fa.fetch(chrom, (start, end))
				fw.write(">{}:{}-{}\n{}\n".format(chrom, start, end, seq))

	elif args.bed_file:
		with open(args.bed_file) as fh:
			for line in fh:
				chrom, start, end = line.strip().split()
				start = int(start) + 1
				end = int(end)
				seq = fa.fetch(chrom, (start, end))
				fw.write(">{}:{}-{}\n{}\n".format(chrom, start, end, seq))

	elif args.regions:
		for region in args.regions:
			chrom, start, end = re.split('[:-]', region)
			start = int(start)
			end = int(end)
			seq = fa[chrom][start-1:end].seq
			fw.write(">{}:{}-{}\n{}\n".format(chrom, start, end, seq))

	else:
		raise Exception("no regions or region file provided")

	if args.out_file:
		fw.close()
	else:
		fw.flush()

def fastx_sample(args):
	fastx_type = fastx_format_check(args.fastx)

	if fastx_type == 'fasta':
		Fastx = pyfastx.Fasta

	elif fastx_type == 'fastq':
		Fastx = pyfastx.Fastq

	else:
		raise Exception("the input file is not fasta or fastq file")

	fx = Fastx(args.fastx)

	if args.num is not None and args.num > 0:
		seq_num = args.num
		if seq_num > len(fx):
			seq_num = len(fx)

	elif args.prop is not None and 0 < args.prop <= 1:
		seq_num = math.ceil(len(fx)*args.prop)

	else:
		raise RuntimeError("specify a right seq number or proportion")

	#set random seed
	if args.seed:
		random.seed(args.seed)
	selected = random.sample(range(len(fx)), k=seq_num)

	if args.out_file is None:
		fw = sys.stdout
	else:
		fw = open(args.out_file, 'w')

	if args.sequential_read:
		writed = 0
		selected = set(selected)
		for it in fx:
			if it.id in selected:
				fw.write(it.raw)
				writed += 1
			else:
				if writed == seq_num:
					break
				_ = it.raw
	else:
		selected.sort()
		for idx in selected:
			fw.write(fx[idx].raw)

	if args.out_file is None:
		fw.flush()
	else:
		fw.close()

def fastx_extract(args):
	fastx_type = fastx_format_check(args.fastx)

	if fastx_type == 'fasta':
		Fastx = pyfastx.Fasta

	elif fastx_type == 'fastq':
		Fastx = pyfastx.Fastq

	else:
		raise Exception("the input file is not fasta or fastq file")

	fx = Fastx(args.fastx)

	if args.out_file:
		fw = open(args.out_file, 'w')
	else:
		fw = sys.stdout

	if args.list_file:
		with open(args.list_file) as fh:
			if args.sequential_read:
				names = {line.strip() for line in fh}
				total = len(names)
				writed = 0

				for it in fx:
					if it.name in names:
						fw.write(it.raw)
						writed += 1
					else:
						if writed == total:
							break
						_ = it.raw
			else:
				for line in fh:
					name = line.strip()
					fw.write(fx[name].raw)

	elif args.names:
		for name in args.names:
			fw.write(fx[name].raw)

	else:
		raise Exception("no sequence name or list file provided")

	if args.out_file:
		fw.close()
	else:
		fw.flush()

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
	parser_build = subparsers.add_parser('index',
		help = "build index for fasta/q file"
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
	parser_info = subparsers.add_parser('stat',
		help = "show detailed statistics information of fasta/q file"
	)
	parser_info.set_defaults(func=fastx_info)
	parser_info.add_argument('fastx',
		help = "fasta or fastq file, gzip support",
		nargs = '+'
	)

	#split command
	parser_split = subparsers.add_parser('split',
		help = "split fasta/q file into multiple files"
	)
	parser_split.set_defaults(func=fastx_split)
	split_group = parser_split.add_mutually_exclusive_group(
		required=True
	)
	split_group.add_argument('-n',
		dest = 'file_num',
		type = int,
		metavar = 'int',
		help = "split a fasta/q file into N new files with even size"
	)
	split_group.add_argument('-c',
		dest = 'seq_count',
		type = int,
		metavar = 'int',
		help = "split a fasta/q file into multiple files containing the same sequence counts"
	)
	parser_split.add_argument('-o', '--out-dir',
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
	parser_fq2fa.add_argument('-o', '--out-file',
		metavar = 'str',
		help = "output file, default: output to stdout"
	)

	parser_fq2fa.add_argument('fastx',
		help = "fastq file, gzip support"
	)

	#get subseq from fasta
	parser_subseq = subparsers.add_parser('subseq',
		help = "get subsequences from fasta file by region"
	)
	parser_subseq.set_defaults(func=fastx_subseq)

	subseq_group = parser_subseq.add_mutually_exclusive_group()
	subseq_group.add_argument('-r', '--region-file',
		metavar = 'str',
		help = "tab-delimited file, one region per line, both start and end position are 1-based"
	)
	subseq_group.add_argument('-b', '--bed-file',
		metavar = 'str',
		help = "tab-delimited BED file, 0-based start position and 1-based end position"
	)
	parser_subseq.add_argument('-o', '--out-file',
		metavar = 'str',
		help = "output file, default: output to stdout"
	)
	parser_subseq.add_argument('fastx',
		help = "input fasta file, gzip support"
	)
	parser_subseq.add_argument('regions',
		help = "format is chr:start-end, start and end position is 1-based, multiple regions were separated by space",
		metavar = 'region',
		nargs = '*'
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
	parser_sample.add_argument('-s', '--seed',
		help = "random seed, default is the current system time",
		type = int,
		default = None,
		metavar = 'int'
	)
	parser_sample.add_argument('--sequential-read',
		action = 'store_true',
		help = "start sequential reading, particularly suitable for sampling large numbers of sequences"
	)
	parser_sample.add_argument('-o', '--out-file',
		metavar = 'str',
		help = "output file, default: output to stdout"
	)
	parser_sample.add_argument('fastx',
		help = "fasta or fastq file, gzip support"
	)

	#extract sequences
	parser_extract = subparsers.add_parser('extract',
		help = "extract full sequences or reads from fasta/q file"
	)
	parser_extract.set_defaults(func=fastx_extract)

	parser_extract.add_argument('-l', '--list-file',
		metavar = 'str',
		help = "a file containing sequence or read names, one name per line"
	)
	parser_extract.add_argument('--reverse-complement',
		help = "output reverse complement sequence",
		action = 'store_true'
	)
	parser_extract.add_argument('--out-fasta',
		help = "output fasta format when extract reads from fastq, default output fastq format",
		action = 'store_true'
	)
	parser_extract.add_argument('-o', '--out-file',
		metavar = 'str',
		help = "output file, default: output to stdout"
	)
	parser_extract.add_argument('--sequential-read',
		action = 'store_true',
		help = "start sequential reading, particularly suitable for extracting large numbers of sequences"
	)
	parser_extract.add_argument('fastx',
		help = "fasta or fastq file, gzip support"
	)
	parser_extract.add_argument('names',
		metavar = 'name',
		help = "sequence name or read name, multiple names were separated by space",
		nargs = '*'
	)

	args = parser.parse_args()

	if hasattr(args, 'func'):
		args.func(args)
	else:
		parser.print_help()

if __name__ == '__main__':
	main()
