import os
import re
import sys
import gzip
import math
import shutil
import random
import pyfastx
import argparse
import traceback
import multiprocessing

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
	for chrom in ids.sort('length', reverse=True):
		idx = min_index(lens)
		fhs[idx].write(fa[chrom].raw)
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

'''
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
'''

def sample_worker(num, fxfile, fxtype, ids, out, lock):
	#create new fasta object in subprocess
	if fxtype == 'fasta':
		fx = pyfastx.Fasta(fxfile)
	elif fxtype == 'fastq':
		fx = pyfastx.Fastq(fxfile)

	if out is None:
		for _id in ids:
			seq = fx[_id].raw
			lock.acquire()
			sys.stdout.write(seq)
			lock.release()
	else:
		with open("{}.{:0>3d}".format(out, num), 'w') as fw:
			for _id in ids:
				fw.write(fx[_id].raw)

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

	if args.threads < 1:
		cpus = 1
	else:
		cpus = args.threads

	selected = random.sample(range(len(fx)), k=seq_num)

	if args.threads > 1:
		#start multiple processes
		count = math.ceil(len(selected)/cpus)
		pool = multiprocessing.Pool(cpus)
		manager = multiprocessing.Manager()
		lock = manager.Lock()

		for i in range(cpus):
			batch = selected[i*count:(i+1)*count]
			pool.apply_async(sample_worker,
				args = (i, args.fastx, fastx_type, batch, args.out_file, lock),
				error_callback = lambda x: print(str(x))
			)

		pool.close()
		pool.join()

		if args.out_file is None:
			sys.stdout.flush()
		else:
			# merge the temp file generated by subprocesses
			with open(args.out_file, 'w') as fw:
				for i in range(cpus):
					temp_file = "{}.{:0>3d}".format(args.out_file, i)

					with open(temp_file) as fh:
						shutil.copyfileobj(fh, fw)

					os.remove(temp_file)

	else:
		if args.out_file is None:
			fw = sys.stdout
		else:
			fw = open(args.out_file, 'w')

		for idx in selected:
			fw.write(fx[idx].raw)

		if args.out_file is None:
			fw.flush()
		else:
			fw.close()


def extract_region_worker(num, fxfile, regions, rc, out, lock):
	fa = pyfastx.Fasta(fxfile)
	rc = '-' if rc else '+'

	if out is None:
		for n, s, e in regions:
			seq = ">{}:{}-{}\n{}\n".format(n, s, e, fa.fetch(n, (s,e), rc))
			lock.acquire()
			sys.stdout.write(seq)
			lock.release()
	else:
		with open("{}.{:0>3d}".format(out, num), 'w') as fw:
			for n, s, e in regions:
				fw.write(">{}:{}-{}\n{}\n".format(n, s, e, fa.fetch(n, (s,e), rc)))

def extract_list_worker(num, fxfile, names, rc, out, lock):
	fa = pyfastx.Fasta(fxfile)

	if out is None:
		if rc:
			for name in names:
				seq = ">{} reverse complement\n{}\n".format(name, fa[name].antisense)
				lock.acquire()
				sys.stdout.write(seq)
				lock.release()
		else:
			for name in names:
				seq = fa[name].raw
				lock.acquire()
				sys.stdout.write(seq)
				lock.release()
	else:
		if rc:
			with open("{}.{:0>3d}".format(out, num), 'w') as fw:
				for name in names:
					fw.write(">{} reverse complement\n{}\n".format(name, fa[name].antisense))
		else:
			with open("{}.{:0>3d}".format(out, num), 'w') as fw:
				for name in names:
					fw.write(fa[name].raw)

def fasta_extract(args):
	fa = pyfastx.Fasta(args.fastx)

	regions = []
	if args.bed_file:
		with open(args.bed_file) as fh:
			for line in fh:
				cols = line.strip().split()
				regions.append((cols[0], int(cols[1])+1, int(cols[2])))
		extract_worker = extract_region_worker

	elif args.region_file:
		with open(args.region_file) as fh:
			for line in fh:
				cols = line.strip().split()
				regions.append((cols[0], int(cols[1]), int(cols[2])))
		extract_worker = extract_region_worker

	elif args.list_file:
		with open(args.list_file) as fh:
			regions = [line.strip().split()[0] for line in fh]
		extract_worker = extract_list_worker

	elif args.region:
		for region in args.region:
			if ':' in region and '-' in region:
				name = region.split(':')[0]
				start, end = region.split(':')[1].split('-')
				regions.append((name, int(start), int(end)))
			else:
				regions.append(region)
		extract_worker = None

	if args.region or args.threads < 1:
		cpus = 1
	else:
		cpus = args.threads

	if cpus == 1:
		if args.out_file:
			fw = open(args.out_file, 'w')
		else:
			fw = sys.stdout

		for region in regions:
			if isinstance(region, tuple):
				n, s, e = region
				if args.reverse_complement:
					seq = fa.fetch(n, (s, e), '-')
				else:
					seq = fa.fetch(n, (s, e))
				fw.write(">{}:{}-{}\n{}\n".format(n, s, e, seq))
			else:
				if args.reverse_complement:
					seq = fa[region].antisense
					fw.write(">{}\n{}\n".format(region, seq))
				else:
					fw.write(fa[region].raw)

		if args.out_file:
			fw.close()
		else:
			fw.flush()

	else:
		count = math.ceil(len(regions)/cpus)
		pool = multiprocessing.Pool(cpus)
		manager = multiprocessing.Manager()
		lock = manager.Lock()

		for i in range(cpus):
			batch = regions[i*count:(i+1)*count]
			pool.apply_async(extract_worker,
				args = (i, args.fastx, batch, args.reverse_complement, args.out_file, lock),
				error_callback = lambda x: print(str(x))
			)

		pool.close()
		pool.join()

		if args.out_file is None:
			sys.stdout.flush()
		else:
			# merge the temp file generated by subprocesses
			with open(args.out_file, 'w') as fw:
				for i in range(cpus):
					temp_file = "{}.{:0>3d}".format(args.out_file, i)

					with open(temp_file) as fh:
						shutil.copyfileobj(fh, fw)

					os.remove(temp_file)

def fastq_worker(num, fxfile, names, outfas, out, lock):
	fq = pyfastx.Fastq(fxfile)

	if out is None:
		if outfas:
			for name in names:
				read = ">{}\n{}\n".format(name, fq[name].seq)
				lock.acquire()
				sys.stdout.write(read)
				lock.release()
		else:
			for name in names:
				read = fq[name].raw
				lock.acquire()
				sys.stdout.write(read)
				lock.release()
	else:
		if outfas:
			with open("{}.{:0>3d}".format(out, num), 'w') as fw:
				for name in names:
					fw.write(">{}\n{}\n".format(name, fq[name].seq))
		else:
			with open("{}.{:0>3d}".format(out, num), 'w') as fw:
				for name in names:
					fw.write(fq[name].raw)

def fastq_extract(args):
	fq = pyfastx.Fastq(args.fastx)

	if args.bed_file:
		with open(args.bed_file) as fh:
			names = [line.strip().split()[0] for line in fh]

	elif args.region_file:
		with open(args.region_file) as fh:
			names = [line.strip().split()[0] for line in fh]

	elif args.list_file:
		with open(args.list_file) as fh:
			names = [line.strip().split()[0] for line in fh]

	elif args.region:
		names = [region for region in args.region]

	if args.region or args.threads < 1:
		cpus = 1
	else:
		cpus = args.threads

	if cpus == 1:
		if args.out_file:
			fw = open(args.out_file, 'w')
		else:
			fw = sys.stdout

		if args.out_fasta:
			if args.out_file is None:
				for name in names:
					sys.stdout.write(">{}\n{}\n".format(name, fq[name].seq))
			else:
				with open(args.out_file, 'w') as fw:
					for name in names:	
						fw.write(">{}\n{}\n".format(name, fq[name].seq))
		else:
			if args.out_file is None:
				for name in names:
					sys.stdout.write(fq[name].raw)
			else:
				with open(args.out_file, 'w') as fw:
					for name in names:
						fw.write(fq[name].raw)

		if args.out_file:
			fw.close()
		else:
			fw.flush()

	else:
		count = math.ceil(len(names)/cpus)
		pool = multiprocessing.Pool(cpus)
		manager = multiprocessing.Manager()
		lock = manager.Lock()

		for i in range(cpus):
			batch = names[i*count:(i+1)*count]
			pool.apply_async(fastq_worker,
				args = (i, args.fastx, batch, args.out_fasta, args.out_file, lock),
				error_callback = lambda x: print(str(x))
			)

		pool.close()
		pool.join()

		if args.out_file is None:
			sys.stdout.flush()
		else:
			# merge the temp file generated by subprocesses
			with open(args.out_file, 'w') as fw:
				for i in range(cpus):
					temp_file = "{}.{:0>3d}".format(args.out_file, i)

					with open(temp_file) as fh:
						shutil.copyfileobj(fh, fw)

					os.remove(temp_file)

def fastx_extract(args):
	if not any((args.region, args.bed_file, args.list_file, args.region_file)):
		raise Exception("No region provided")

	fastx_type = fastx_format_check(args.fastx)

	if fastx_type == 'fasta':
		fasta_extract(args)

	elif fastx_type == 'fastq':
		fastq_extract(args)

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

	'''
	#get subseq from fasta
	parser_subseq = subparsers.add_parser('subseq',
		help = "get subsequences from fasta file by id or name with region"
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
	'''

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
	parser_sample.add_argument('-o', '--out-file',
		metavar = 'str',
		help = "output file, default: output to stdout"
	)
	parser_sample.add_argument('-t', '--threads',
		metavar = 'int',
		type = int,
		default = 1,
		help = "number of CPUs (or processes) to launch, default: 1"
	)
	parser_sample.add_argument('fastx',
		help = "fasta or fastq file, gzip support"
	)

	#extract sequences
	parser_extract = subparsers.add_parser('extract',
		help = "extract subsequences or reads from fasta/q file"
	)
	parser_extract.set_defaults(func=fastx_extract)
	
	extract_group = parser_extract.add_mutually_exclusive_group()
	extract_group.add_argument('-l', '--list-file',
		metavar = 'str',
		help = "a file containing sequence or read names, one name per line"
	)
	extract_group.add_argument('-b', '--bed-file',
		metavar = 'str',
		help = "tab-delimited BED file, 0-based start position and 1-based end position"
	)
	extract_group.add_argument('-r', '--region-file',
		metavar = 'str',
		help = "tab-delimited file, but both start and end position are 1-based"
	)
	'''
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
	'''
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
	parser_extract.add_argument('-t', '--threads',
		metavar = 'int',
		type = int,
		default = 1,
		help = "number of CPUs (or processes) to launch, default: 1"
	)
	parser_extract.add_argument('fastx',
		help = "fasta or fastq file, gzip support"
	)
	parser_extract.add_argument('region',
		help = "format is chr or chr:start-end, multiple regions were separated by space",
		nargs = '*'
	)

	args = parser.parse_args()

	if hasattr(args, 'func'):
		args.func(args)
	else:
		parser.print_help()

if __name__ == '__main__':
	main()
