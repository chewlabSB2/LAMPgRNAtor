#!/usr/bin/env python3

import os
import argparse
from LAMPgRNAtor.LAMP_utils import *
import GLAPD
from LAMPgRNAtor._version import __version__ as version

DESCRIPTION = '''LAMP Primer Designer'''

def get_args():
	parser = argparse.ArgumentParser(
		prog = f'{version}',
		description = DESCRIPTION
	)
	parser.add_argument('-r', '--reference', dest = 'reference', required = True, metavar="fasta", 
						help="Reference for Guide Score Calculation")
	parser.add_argument('-p', '--prefix', default = 'LAMP_primers', 
						help='file name prefix to your file names (Default: %(default)s)')	
	
	args = parser.parse_args()

	return args

def main():
	args = get_args()
	ref = get_sequence(args.reference)
	GLAPD_prefix = args.prefix + '.GLAPD'
	GLAPD.single_candidate(str.encode(ref), str.encode(args.prefix), STORE_PATH)
	GLAPD.design_LAMP(str.encode(args.prefix), str.encode(ref), str.encode(GLAPD_prefix), STORE_PATH)

if __name__ == '__main__':
	main()
