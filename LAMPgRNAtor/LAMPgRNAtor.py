#!/usr/bin/env python3

import os
import sys
import logging
import logging.handlers
import time
import argparse
import tqdm
from LAMPgRNAtor.LAMP_utils import *
import GLAPD
from LAMPgRNAtor._version import __version__ as version

logger = logging.getLogger(__name__)
#version = "LAMPgRNAtor v1.0"

DESCRIPTION = f'''
{version} combine module Predicts for Cas12a gRNAs \n
{version} scores for Conservation and Entropy \n 
{version} Utilizes CPU cores under the threads option \n
'''

def get_args():
	parser = argparse.ArgumentParser(
		prog = f'{version}',
		description = DESCRIPTION
	)
	parser.add_argument('-r', '--reference', dest = 'reference', metavar="fasta", 
						help="Reference for Guide Score Calculation")
	parser.add_argument('-m','--MSA', dest = "MSA", metavar="fasta", 
						help="Sequences to align")
	parser.add_argument('-bc', '--bowtie-conservation', dest = "bowtie", 
						help="Path to fasta file for Bowtie conservation calculation")
	parser.add_argument('-bo', '--bowtie-offtarget', dest = "offtarget", 
						help="Path to fasta file for Bowtie offtarget  Filtering")
	
	parser.add_argument('-p', '--prefix', default = 'LAMPgRNAtor', 
						help='file name prefix to your file names (Default: %(default)s)')	
	parser.add_argument('-d', '--debug',
						help='Print lots of debugging statements',
						action="store_const",dest="loglevel",const=logging.DEBUG,
						default=logging.INFO)
	parser.add_argument('-c', '--concatenate', default=None, dest = 'concat',
						help='Concatenate between F1c and B1c eg. TTTT [make sure it is A, C, G, T] (Default: %(default)s)')
	parser.add_argument('-l', '--LAMP',  dest='LAMP_csv', nargs = '+',
						metavar="txt", help="accepts a list of GLAPD and LAMP output")
	parser.add_argument('-t', '--threshold', default=50, help="Minimum cut-off point for DeepCpf1 score (Default: %(default)s)")

	parser.add_argument('--no-GLAPD', dest='find_primers', action="store_false", 
						help='Find LAMP Primers using GLAPD (Default: True)')
	parser.add_argument('--overlap', action="store_false", 
						help='Allow overlap of gRNA with F1c and B1c regions (Default: (No Overlap))')
	parser.add_argument('--threads', default=4, type=int, dest = 'CPU', 
						help="Number of Core to use (Default: %(default)s)")
	parser.add_argument('--mismatch', default=5, type=int,  choices=range(0,10),
						help='Number of mismatches allowed for Consensus Search (Default: %(default)s)')
	parser.add_argument('--conservation-index', dest="conservation_index", default=None, 
						help='Input an index for conservation calculation if you have!')
	parser.add_argument('--offtarget-index', dest="offtarget_index", default=None, 
						help='Input an index for offtargets if you have!')	
	parser.add_argument('--keep-tmp', action="store_false", dest='keep_temp', 
						help='Keep Temporary files (Default: Does not Keep)')
	parser.add_argument('--no-plot', action="store_false",dest="plot", 
						help='Does not Plot graph')
	parser.add_argument('--no-gRNA', action='store_false', dest='find_gRNA', help='Switch off gRNA search mode' )
	parser.add_argument('--ignore-loop', action='store_false', dest='LOOP', help='ignore LOOPs, doesn\' matter if LF or LB exists' )
	
	args = parser.parse_args()
	if args.offtarget_index or args.conservation_index:
		args.keep_temp = False

	return args

def calculate_MSA(args):
	consensus, conservation, entropy = None, None, None
	MSA_ed, do_bowtie = False, False
	if args.MSA: 
		aln_array, MSA = read_alignment(True, args.MSA)
		logger.debug(f"There are {aln_array} aligned sequences")
		if MSA:
			if len(aln_array) < 100: 
				logger.debug("Number of Sequences insufficient!!!")
				logger.debug("Conservation and Entropy Score might not be accurate")
			logger.info("Will Calculate Conservation and Entropy Score")
			consensus, conservation, entropy = main_calc(aln_array, args.prefix, plot = args.plot)
			MSA_ed = True
		else:
			logger.info("Will Bowtie to get Conservation Score") 
			do_bowtie = True
	
	if args.bowtie:
		do_bowtie = True

	return consensus, (conservation, entropy), (MSA_ed, do_bowtie)

def bowtie_main(args, crRNA_fasta, mode = 'offtarget'):
	if mode == 'offtarget' and not args.offtarget:
		return {}, []

	temp_files_to_remove = []
	bowtie_index_prefix = None
	bowtie_path, bowtie_build_path = check_bowtie(args)
	
	if mode == 'offtarget':
		_prefix = args.prefix  + '_offtarget' 
		bowtie_reference = args.offtarget
		samfile = args.prefix + '_offtarget.sam'
		bowtie_index_prefix = args.offtarget_index if args.offtarget_index else None
	elif mode == 'conservation':
		_prefix = args.prefix + '_conservation'
		bowtie_reference = args.bowtie
		samfile = args.prefix + '_conservation.sam'
		bowtie_index_prefix = args.conservation_index if args.conservation_index else None

	if not bowtie_index_prefix:
		logger.info(f"Building Bowtie Index for {mode}! This might take a while")
		bowtie_build_cmdline, bowtie_index_prefix = bowtie_build_cmd(bowtie_reference, _prefix, bowtie_build_path, args.CPU)
		success_bowtie_build = run_shell_command(bowtie_build_cmdline)
		if not success_bowtie_build:
			raise AlignmentError("Error during bowtie-build")
		temp_files_to_remove.append(_prefix + '_output_build.txt')
		temp_files_to_remove.append(_prefix + '_error_build.txt')

	#crRNA_fasta == Reference Fasta File
	logger.info(f"Bowtie Starting for {mode}!")
	mismatch = args.mismatch if args.mismatch <= 3 else 3
	cmdline_bowtie = bowtie_cmd(bowtie_index_prefix, crRNA_fasta, _prefix, mismatch, bowtie_path, args.CPU)
	success_bowtie = run_shell_command(cmdline_bowtie)
	if not success_bowtie:
		raise AlignmentError(f"Error during {mode} Scoring Alignment")

	bowtie_summary = bowtie_to_summary(samfile, mode = mode)
	logger.info(f"Bowtie Done for {mode}!")
	temp_files_to_remove.append(samfile)

	onlyfiles = [f for f in os.listdir() if f.startswith(f"{_prefix}_alignment.index")]
	temp_files_to_remove += onlyfiles

	return bowtie_summary, temp_files_to_remove

def _create_fasta(primer_class, prefix):
	LAMP_handle = prefix + '_LAMP.fasta'
	with open(LAMP_handle,'w') as f:
		for LAMP in primer_class:
			for p in PRIMER_LIST:
				primer = getattr(LAMP, p)
				if not primer: continue
				primer_seq = primer.seq
				if primer.reverse:
					primer_seq = reverseC(primer_seq)
				f.write(f'>{primer.id}\n{primer_seq}\n')

			for g in LAMP.gRNAs:
				gRNA_seq = g.seq
				if g.reverse:
					gRNA_seq = reverseC(gRNA_seq)
				f.write(f'>{g.id}\n{gRNA_seq}\n')

	return LAMP_handle

def _update_LAMP(bowtie_summary, LAMP_Class, mode = 'conservation'):
	for LAMP in LAMP_Class:
		for primer in PRIMER_LIST:
			current_primer = getattr(LAMP, primer)
			if not current_primer: continue
			current_primer_id = current_primer.id
			if current_primer_id in list(bowtie_summary.keys()):
				if mode == 'conservation':
					current_primer.c_score_bowtie = bowtie_summary[current_primer_id]
				elif mode == 'offtarget':
					for o in bowtie_summary[current_primer_id]:
						current_primer.offtargets.append(o)

		for g in LAMP.gRNAs:
			gRNA_id = g.id
			if gRNA_id in list(bowtie_summary.keys()):
				if mode == 'conservation':
					g.c_score_bowtie = bowtie_summary[gRNA_id]
				elif mode == 'offtarget':
					for o in bowtie_summary[gRNA_id]:
						g.offtargets.append(o)

def _write_LAMP(args, LAMP_Class, status):
	bowtie_ed, MSA_ed, offtarget_ed = status
	if not os.path.isdir(args.prefix + RESULT_DIR): 
		os.mkdir(args.prefix + RESULT_DIR)  

	os.chdir(args.prefix + RESULT_DIR)
	for LAMP in LAMP_Class:
		LAMP_id = LAMP.id
		with open(f'{args.prefix}_{LAMP_id}.csv', 'w') as write_all:
			
			## Write Finalize Sequence (DONE)
			LAMP.create_inner_primers(args.concat)
			write_all.write(f'{LAMP_id} found between {LAMP.start} and {LAMP.end}\n')
			for primer in PRIMER_LIST_WRITE:
				p = getattr(LAMP, primer)
				if not p: continue
				seq = p if primer in INNER_PRIMERS else p.seq
				write_all.write(f'{primer},{seq}\n')
			for i,g in enumerate(LAMP.gRNAs):
				reverse_state = 'rev' if g.reverse else 'fow'
				write_all.write(f'gRNA-{i}-{reverse_state},{g.seq}\n')
			write_all.write('\n\n')

			## Write Scores
			write_all.write('MAIN SUMMARY')
			write_all.write(",".join(HEADER_MAIN_1))
			if MSA_ed:
				write_all.write("," + ",".join(HEADER_MSA))
			if bowtie_ed:
				write_all.write("," + HEADER_BOWTIE_CONSERVATION)
			if offtarget_ed:
				write_all.write("," + HEADER_BOWTIE_OFFTARGETS)

			write_all.write('\n')
			for primer in PRIMER_LIST:
				p = getattr(LAMP, primer)
				if not p: continue
				sequence = reverseC(p.seq) if p.reverse else p.seq
				write_all.write(f'{p.id},{sequence},{p.pos},{p.len},{p.reverse}')
				if MSA_ed:
					write_all.write(f',{p.found},{p.consensus_query},{p.mismatch},{p.cigar}')
					write_all.write(f',{p.e_sum},{p.c_mean},{p.c_SD},0')
				if bowtie_ed:
					write_all.write(f',{p.c_score_bowtie}')
				if offtarget_ed and p.offtargets:
					for genome in p.offtargets:
						write_all.write(f",{genome}")
				write_all.write('\n')

			if LAMP.gRNAs:
				for g in LAMP.gRNAs:
					write_all.write(f'{g.id},{g.seq},{g.pos},{g.len},{g.reverse}')
					if MSA_ed:
						write_all.write(f',{g.found},{g.consensus_query},{g.mismatch},{g.cigar}')
						write_all.write(f',{g.e_sum},{g.c_mean},{g.c_SD},{g.score}')
					if bowtie_ed:
						write_all.write(f',{g.c_score_bowtie}')
					if offtarget_ed and g.offtargets:
						for genome in g.offtargets:
							write_all.write(f",{genome}")
					write_all.write('\n')

			write_all.write('\n\n')
			## Write Conservation Supplementary Score
			if MSA_ed:
				write_all.write('CONSERVATION SCORES\n')
				for primer in PRIMER_LIST:
					p = getattr(LAMP, primer)
					if not p: continue
					if not p.c_score_list: continue
					conservation_list = [str(round(x, 5)) for x in p.c_score_list]
					write_all.write(f'{p.id},' + ",".join(conservation_list) + "\n")

				if LAMP.gRNAs:
					for g in LAMP.gRNAs:
						if not g.c_score_list: continue
						conservation_list = [str(round(x, 5)) for x in g.c_score_list]
						write_all.write(f'{g.id},' + ",".join(conservation_list) + "\n")

				write_all.write('\n\n')
				write_all.write('ENTROPY SCORES\n')
			## Write Entropy Supplementary Score
				for primer in PRIMER_LIST:
					p = getattr(LAMP, primer)
					if not p: continue
					if not p.e_score_list: continue
					entropy_list = [str(round(x, 5)) for x in p.e_score_list]
					write_all.write(f'{p.id},' + ",".join(entropy_list) + "\n")

				if LAMP.gRNAs:
					for g in LAMP.gRNAs:
						if not g.e_score_list: continue
						conservation_list = [str(round(x, 5)) for x in g.e_score_list]
						write_all.write(f'{g.id},' + ",".join(entropy_list) + "\n")

				write_all.write('\n\n')

			write_all.write('ALTERNATIVES\n')
			## Write Alternatives
			write_all.write(",".join(HEADER_MAIN_1))
			write_all.write("," + ",".join(ALTERNATIVE_MAIN) + '\n')
			for primer in PRIMER_LIST:
				p = getattr(LAMP, primer)
				if not p: continue
				if not p.alternative: continue
				for a in p.alternative:
					write_all.write(f'{p.id},{reverseC(p.seq) if p.reverse else p.seq},{p.pos},{p.len},{p.reverse}')
					write_all.write(f',{a.query},{a.pos[0]},{a.editDistance},{a.cigar}\n')

			for g in LAMP.gRNAs:
				if g.alternative:
					for a in g.alternative:
						write_all.write(f'{g.id},{reverseC(g.seq) if g.reverse else g.seq},{g.pos},{g.len},{g.reverse}')
						write_all.write(f',{a.query},{a.pos[0]},{a.editDistance},{a.cigar}\n')

	os.chdir("..")

def debugging(primer_class, message = "Test"):
	## FOR DEBUGGING PURPOSES
	logger.info(f"Debug in Stage: {message}")
	for p in primer_class:
		logger.info(f"{p.id} F3: {p.F3.seq}")
		logger.info(f"{p.id} F2: {p.F2.seq}")
		logger.info(f"{p.id} F1c: {p.F1c.seq}")
		logger.info(f"{p.id} F1c: {p.F1c.pos}")
		logger.info(f"{p.id} F1c: {p.F1c.len}")
		logger.info(f"{p.id} B3: {p.B3.seq}")
		logger.info(f"{p.id} B2: {p.B2.seq}")
		logger.info(f"{p.id} B1c: {p.B1c.seq}")
		logger.info(f"{p.id} B1c: {p.B1c.pos}")
		logger.info(f"{p.id} B1c: {p.B1c.len}")
		if p.LF: logger.info(f"{p.id} LF: {p.LF.seq}")
		if p.LB: logger.info(f"{p.id} LB: {p.LB.seq}")

def main():
	args = get_args()
	start_time = time.time()
	logging.basicConfig(level=args.loglevel, format=FORMAT)
	logger.info(f"LAMPgRNAtor v{version} starting!!")
	temp_files_to_remove = []
	bowtie_ed = MSA_ed = offtarget_status = False

	## Reference for 1) to find gRNA 2) LAVA
	reference = get_sequence(args.reference)
	reference = Wavelet_Tree(reference)
	
	primer_class = []
	if args.find_primers:
		store_path = str.encode(os.getcwd() + '/')
		logger.info(f"Finding Potential Primers!")
		minLen = 18 if args.find_gRNA else 0
		GLAPD_prefix = args.prefix + '.GLAPD'
		ref = reference.ReconstructSequence()
		GLAPD.single_candidate(str.encode(ref), str.encode(args.prefix), store_path)
		
		logger.info(f"Finding LAMP Primer Sets!")
		GLAPD.design_LAMP(str.encode(args.prefix), str.encode(ref), str.encode(GLAPD_prefix), store_path, minlen = minLen)
		
		del ref
		args_stored = (args.overlap, args.threshold, [GLAPD_prefix])
		primer_class = open_LAMP_files(args_stored, reference, primer_class)

	## Open LAMP Primer Files 
	if args.LAMP_csv:
		logger.info(f"Opening LAMP Files!")
		args_stored = (args.overlap, args.threshold, args.LAMP_csv)
		primer_class = open_LAMP_files(args_stored, reference, primer_class)
		#debugging(primer_class)

	logger.info(f"Total Number of LAMP Primer Sets: {len(primer_class)}")
	
	## Calculate Conservation & Entropy
	logger.info(f"Calculating MSA, might take awhile!")
	consensus, score, MSA_info = calculate_MSA(args)
	MSA_ed, do_bowtie = MSA_info
	if MSA_ed: 
		consensus = Wavelet_Tree(consensus)

	## Find gRNAs and/or Calculate Scores for each Primer
	logger.info(f"Searching for Primers in Consensus Sequence!")
	primer_class = multi_search(args, consensus, reference, primer_class, MSA_ed = MSA_ed, scores = score)
	del score
	del consensus
	del reference

	#debugging(primer_class, message = "Updated")
	## Create Fasta for Bowtie 
	LAMP_handle = _create_fasta(primer_class, args.prefix)
	conservation_summary, offtarget_summary = None, None

	if sys.platform != 'win32':
		## Conservation Bowtie
		logger.info(f"Bowtie for Conservation!")
		conservation_summary, temp_files_to_append = bowtie_main(args, LAMP_handle, mode = 'conservation')
		temp_files_to_remove += temp_files_to_append

		## Offtarget Bowtie
		logger.info(f"Bowtie for Offtarget!")
		offtarget_summary, temp_files_to_append = bowtie_main(args, LAMP_handle, mode = 'offtarget')
		temp_files_to_remove += temp_files_to_append
	
	else:
		''' FUTURE UPDATE '''
		logger.info(f"Skipping Bowtie!")
	
	## Update LAMP Objects - Offtarget & Conservation
	if conservation_summary: 
		_update_LAMP(conservation_summary, primer_class, mode = 'conservation')
		bowtie_ed = True
	
	if offtarget_summary: 
		_update_LAMP(offtarget_summary, primer_class, mode = 'offtarget')
		offtarget_status = True

	## Write Primers
	logger.info(f"Writing Output!")
	status = (bowtie_ed, MSA_ed, offtarget_status)
	_write_LAMP(args, primer_class, status)

	if args.keep_temp:
		for fname in temp_files_to_remove:
			os.remove(fname)      
		logger.info("Temporary files removed!")

	end_time = time.time()
	total_time = round(end_time - start_time,3)
	logger.info(f"{version} has Ended in {total_time} seconds!")

if __name__ == '__main__':
	main()