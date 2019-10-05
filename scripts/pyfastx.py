import pyfastx
import argparse

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
		help = "Split fasta "
	)


