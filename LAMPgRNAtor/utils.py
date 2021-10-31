#!/usr/bin/env python3

import os
import subprocess
import sys
from textwrap import dedent
import shlex
import logging
import logging.handlers
from Bio import AlignIO, SeqIO
import numpy as np
import matplotlib.pyplot as plt



FORMAT = '%(asctime)s - %(module)-16s - %(levelname)s - %(message)s'
PYTHON_FILE = os.path.dirname(os.path.abspath(__file__))
shquote = shlex.quote

ALPHABET = ['A', 'C', 'G', 'T', 'N']

class ExecutionError(Exception):
	pass

class AlignmentError(Exception):
	pass

def reverseC(seq):
	RT = {'A':'T','C':'G','G':'C','T':'A', 'N':'N'}
	reverseComplement = ''
	for i in seq:
		nt = RT.get(i)
		reverseComplement += nt
	return reverseComplement[::-1]

def get_sequence(reference):
	record = SeqIO.read(reference, "fasta")
	return str(record.seq).upper()

##-------------------------------------------------------------------------------------------------
## Bowtie Offtarget and Conservation Scoring

def is_tool(name):
	## Check whether `name` is on PATH and marked as executable.
	from shutil import which
	return which(name) is not None

def check_bowtie(args):
	if not is_tool('bowtie') or not is_tool('bowtie-build'):
		if not is_tool(PYTHON_FILE + '/misc/bowtie') or not is_tool(PYTHON_FILE + '/misc/bowtie-build'): 
			raise AlignmentError("Bowtie is not installed")

	if args.bowtie or args.bowtie_index or args.gRNA_mode[0].lower() == 'low': #Edit
		if not is_tool('bowtie') or not is_tool('bowtie-build'):
			if not is_tool(PYTHON_FILE + '/misc/bowtie') or not is_tool(PYTHON_FILE + '/misc/bowtie-build'): 
				raise AlignmentError("Bowtie is not installed")
			else:
				bowtie_path = PYTHON_FILE + '/misc/bowtie'
				bowtie_build_path = PYTHON_FILE + '/misc/bowtie-build'
		else:           
			bowtie_path = 'bowtie'
			bowtie_build_path = 'bowtie-build'

	return bowtie_path, bowtie_build_path

		
def bowtie_build_cmd(bowtie_reference, prefix, bowtie_path = 'bowtie-build', thread = 4):
	cmd = '%s --threads %d %s %s_alignment.index 1>>%s_output_build.txt 2>>%s_error_build.txt'%(bowtie_path, thread, bowtie_reference, prefix, prefix, prefix)
	return cmd, prefix + '_alignment.index'

def bowtie_cmd(bow, reference_fasta, prefix, mismatch, bowtie_path = 'bowtie', thread = 4):
	cmd = '%s -a -S %s -f %s %s.sam -v %d --threads %d'%(bowtie_path, bow, reference_fasta, prefix, mismatch, thread)
	return cmd

##-------------------------------------------------------------------------------------------------
##Conservation And Entropy Score Calculation

def read_alignment(MSA, handle):
	## Return MSA if alignment works else return Null
	if MSA:
		try:
			return AlignIO.read(handle, 'fasta'), True
		except:
			try:
				test = SeqIO.parse(handle, 'fasta')
				return None , False
			except Exception as error:
				raise AlignmentError("\nERROR: Problem reading in {}: {}".format(handle, str(error)))
	else:
		try:
			test = AlignIO.read(handle, 'fasta')
			return None, True
		except:
			try:
				test = SeqIO.parse(handle, 'fasta')
				return None, False
			except Exception as error:
				raise AlignmentError("\nERROR: Problem reading in {}: {}".format(handle, str(error)))

def prettify_gaps(MSA):
	## Converts all bases to uppercase
	## Replace all gaps by 'N' in all sequences in the alignment post-alignment	
	from Bio.Seq import Seq
	for seq_ in MSA:
		_seq = str(seq_.seq)
		_seq = _seq.upper()
		_seq = _seq.replace('-', 'N')
		seq_.seq = Seq(_seq)

def calc_frequencies(MSA):
	## Returns a matrix of nucleotide frequecy for each position for entropy and conservation score
	aln_array = np.array(MSA) 
	mat = np.zeros((len(ALPHABET), aln_array.shape[1]))
	for ni, nuc in enumerate(ALPHABET):
		mat[ni] = np.mean(aln_array == nuc, axis=0)
	mat /= np.sum(mat, axis=0) 
	return mat

def calc_entropy(mat):
	## Shannon's Entropy for each position in MSA
	return [round(x, 5) for x in (-1) * np.sum(mat * np.log2(mat + 1e-15), axis=0)]

def consensus_seq(mat):
	## Get Consensus Sequence from matrix for wessels gRNA scoring
	consensus = []
	## max value from each array
	consensus_value = np.max(mat, axis=0) 
	## position/index of max value in array
	a = np.argmax(mat, axis=0) 
	
	for i in a:
		consensus.append(ALPHABET[i])

	return ''.join(consensus), consensus_value

def convert_to_fasta(consensus_sequence, prefix):
	from Bio.Seq import Seq
	from Bio.SeqRecord import SeqRecord

	sequence_to_save = SeqRecord(seq = Seq(consensus_sequence), id = 'Consensus', description = 'Consensus_Sequence')
	handle = open(prefix + "_consensus_sequence.fasta","w")
	SeqIO.write(sequence_to_save,handle,"fasta")
	handle.close()

def write_score(entropy_score_table, conservation_score_table, consensus_sequence, prefix):
	csv_file = prefix + '_position_score.csv'
	score_csv = open(csv_file, 'w')
	score_csv.write('Position\tEntropy\tConservation\tBase Pair\n')
	for i in range(len(entropy_score_table)):
		score_csv.write(str(i+1) + '\t' + str(entropy_score_table[i]) + '\t' + str(conservation_score_table[i]) +  '\t' + consensus_sequence[i] + '\n')
	score_csv.close()

def get_moving_average(score, window_size = 10):
	if window_size == 0:
		return score
	
	i = 0
	moving_averages = []
	if window_size == 0: window_size = 1
	while i < len(score) - window_size + 1:
		this_window = score[i : i + window_size]
		window_average = sum(this_window) / window_size
		moving_averages.append(window_average)
		i += 1
	
	return moving_averages

def plot_conservation(entropy, conservation, prefix):
	fig, ax = plt.subplots(figsize=(15, 5), sharex=True)

	ax.plot([i for i in range(len(conservation))], conservation, linestyle="", color="blue", marker="o", markersize = 0.5)
	ax.set_ylabel("Conservation %",color="blue",fontsize=14)
	ax.tick_params(axis="x", labelsize=12) 
	ax.tick_params(axis='x', which='both', labelsize=10, labelbottom=True)
	
	ax2 = ax.twinx()
	ax2.plot([i for i in range(len(entropy))], entropy, linestyle="", color="red", marker="o", markersize = 0.5)
	ax2.set_ylabel("Entropy",color="red",fontsize=14)

	fig.suptitle(f'{prefix} Entropy | Conservation Score',  fontsize=16)
	figure = plt.gcf() # get current figure
	figure.set_size_inches(12, 8)
	plt.grid(False)
	plt.savefig(prefix + '_combined_genome_plot.png', dpi = 360)

def main_calc(MSA, prefix, plot = True):
	prettify_gaps(MSA)
	matrix = calc_frequencies(MSA)
	entropy_list = calc_entropy(matrix)
	consensus, conservation_list = consensus_seq(matrix)
	conservation_list = [round(i, 5) for i in conservation_list]
	write_score(entropy_list, conservation_list, consensus, prefix)
	convert_to_fasta(consensus, prefix)
	if plot:
		for ma in [0, 10, 25, 50, 100, 250]:
			adjusted_entropy = get_moving_average(entropy_list, ma)
			adjusted_conservation = get_moving_average(conservation_list, ma)
			plot_conservation(adjusted_entropy, adjusted_conservation, f'{prefix}-{ma}')
	return consensus, conservation_list, entropy_list 

def plot_conservation(entropy, conservation, prefix):
	fig, ax = plt.subplots(figsize=(15, 5), sharex=True)

	ax.plot([i for i in range(len(conservation))], conservation, linestyle="", color="blue", marker="o", markersize = 0.5)
	ax.set_ylabel("Conservation %",color="blue",fontsize=14)
	ax.tick_params(axis="x", labelsize=12) 
	ax.tick_params(axis='x', which='both', labelsize=10, labelbottom=True)
	
	ax2 = ax.twinx()
	ax2.plot([i for i in range(len(entropy))], entropy, linestyle="", color="red", marker="o", markersize = 0.5)
	ax2.set_ylabel("Entropy",color="red",fontsize=14)

	fig.suptitle(f'{prefix} Entropy | Conservation Score',  fontsize=16)
	figure = plt.gcf() # get current figure
	figure.set_size_inches(12, 8)
	plt.grid(False)
	plt.savefig(prefix + '_combined_genome_plot.png', dpi = 360)

##-------------------------------------------------------------------------------------------------
##Other Functions

def print_error(message, **kwargs):
	"""
	Formats *message* with *kwargs* using :meth:`str.format` and
	:func:`textwrap.dedent` and uses it to print an error message to
	``sys.stderr``.
	"""
	print("\nERROR: " + dedent(message.format(**kwargs)).lstrip("\n")+"\n", file = sys.stderr) 

def run_shell_command(cmd, raise_errors = False, extra_env = None):
	"""
	Run the given command string via Bash with error checking.
	Returns True if the command exits normally.  Returns False if the command
	exits with failure and "raise_errors" is False (the default).  When
	"raise_errors" is True, exceptions are rethrown.

	If an *extra_env* mapping is passed, the provided keys and values are
	overlayed onto the default subprocess environment.
	"""
	env = os.environ.copy()

	if extra_env:
		env.update(extra_env)

	shargs = ['-c', "set -euo pipefail; " + cmd]

	if os.name == 'posix':
		shellexec = ['/bin/bash']
	else:
		# We try best effort on other systems. For now that means nt/java.
		shellexec = ['env', 'bash']

	try:
		# Use check_call() instead of run() since the latter was added only in Python 3.5.
		subprocess.check_output(
			shellexec + shargs,
			shell = False,
			stderr = subprocess.STDOUT,
			env = env)

	except subprocess.CalledProcessError as error:
		print_error(
			"{out}\nshell exited {rc} when running: {cmd}{extra}",
			out = error.output,
			rc  = error.returncode,
			cmd = cmd,
			extra = "\nAre you sure this program is installed?" if error.returncode==127 else "",
		)
		if raise_errors:
			raise
		else:
			return False

	except FileNotFoundError as error:
		print_error(
			"""
			Unable to run shell commands using {shell}!
			Module requires {shell} to be installed.
			""",
			shell = ' and '.join(shellexec)
		)
		if raise_errors:
			raise
		else:
			return False

	else:
		return True  


