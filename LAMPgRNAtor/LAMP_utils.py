#!/usr/bin/env python3

## LAMP Class
import logging
import logging.handlers
import re
import os
from collections import defaultdict
from matplotlib.lines import Line2D
from LAMPgRNAtor.utils import * 
from LAMPgRNAtor.simpleSearch import *

logging.getLogger('matplotlib').setLevel(logging.WARNING)

logger = logging.getLogger(__name__)
PRIMER_LIST = ['F3', 'B3', 'F2', 'B2', 'F1c', 'B1c', 'LF', 'LB']
PRIMER_LIST_REVERSE = ['B3', 'B2', 'F1c', 'LF']
THREE = ['F3', 'B3']
TWO = ['F2','B2']
ONE = ['F1c', 'B1c']
LOOP = ['LB', 'LF']
LAVA_LOOP = ['FLOOP', 'BLOOP']
PRIMER_LIST_WRITE = ['F3', 'B3', 'FIP', 'BIP', 'LF', 'LB']
INNER_PRIMERS = ['FIP', 'BIP']

PLOT_DIR = '_LAMPgRNAtor_plots'
RESULT_DIR = '_LAMPgRNAtor_csv'

HEADER_MAIN_1 = [
	'Name',
	'Sequence',
	'Pos',
	'Length',
	'Reverse',
]

HEADER_MSA = [
	'Found',
	'Sequence_Found',
	'Mismatch',
	'Cigar',
	'Entropy_Sum',
	'Conservation_Mean_Score',
	'Conservation_SD_Score',
]

ALTERNATIVE_MAIN = [
	'Consensus_Query',
	'MatchPos',
	'Mismatch',
	'Cigar',
]

HEADER_BOWTIE_OFFTARGETS = 'offtargets'
HEADER_BOWTIE_CONSERVATION = 'Conservation%'

STORE_PATH = str.encode(os.getcwd() + '/')

WINDOW_SIZE = 25
gRNA_LENGTH = 27

##-------------------------------------------------------------------------------------------------
## LAMP Class & Internal Functions

class LAMP_set():

	def __init__(self, LAMP_id, overlap = False, threshold = 0):
		self.id = LAMP_id
		self.F3 = None
		self.B3 = None
		self.F2 = None
		self.B2 = None
		self.F1c = None
		self.B1c = None
		self.LF = None
		self.LB = None

		## Inner Primers
		self.BIP = None
		self.FIP = None

		## gRNAs
		self.start_gRNA = 0
		self.end_gRNA = float('inf')
		self.gRNAs = []
		self.overlap = overlap
		self.sequence_gRNA = None
		self.short = False
		self.threshold = threshold

		self.start = -1 ## F3 start
		self.end = float('inf') ## B3 end
		self.c_score_list = []
		self.e_score_list = []

	def create_inner_primers(self, concat = None):
		if concat:
			self.FIP = self.F1c.seq + concat + self.F2.seq
			self.BIP = self.B1c.seq + concat + self.B2.seq
		else:
			self.FIP = self.F1c.seq + self.F2.seq
			self.BIP = self.B1c.seq + self.B2.seq

	def _update_start_end(self, max_length = 400):
		self.start = self.F3.consensus_pos[0]
		self.end = self.B3.consensus_pos[1] #+ self.B3.len

		## ADD a Warning check 

	def get_Cas12gRNAs(self, ref, deepCpf1):
		self._get_gRNA_sequence(ref)
		self._generate_gRNAs()
		self._predict_gRNAs(deepCpf1)

	def _get_gRNA_sequence(self, ref):
		self.start_gRNA = self.F1c.pos if self.overlap else (self.F1c.pos + self.F1c.len)
		self.end_gRNA = (self.B1c.pos + self.B1c.len) if self.overlap else self.B1c.pos
		if self.start_gRNA != None and self.end_gRNA != None:
			if self.start_gRNA > len(ref) or self.end_gRNA > len(ref):
				raise ExecutionError('Sequence provided is too short')

		self.sequence_gRNA = ref[self.start_gRNA: self.end_gRNA]
		if len(self.sequence_gRNA) < 27:
			self.short = True #Skip Search
	   
	def _generate_gRNAs(self):
		if self.short:
			return

		count = 0
		PAM_nucleotide = ['A', 'G', 'C']
		pam = []
		for a in PAM_nucleotide:
			pattern = 'TTT' + a
			pam = bitapSearch(self.sequence_gRNA, pattern, mismatch = 0)
			if pam: 
				for p in pam:
					if ((p - 4) < 0) or ((p+30)>len(self.sequence_gRNA)): continue
					pos = self.start_gRNA + p
					seq = self.sequence_gRNA[p-4:p+30]
					if len(seq) != 34: continue
					count += 1
					gRNA_temp = Primer(f'gRNA', f'{self.id}-gRNA_{count}', seq, pos, source = 'gRNA', reverse = False)
					self.gRNAs.append(gRNA_temp)

		PAM_nucleotide = ['G', 'C', 'T']
		pam = []
		for a in PAM_nucleotide:
			pattern = a + 'AAA'
			pam = bitapSearch(self.sequence_gRNA, pattern, reverse = True, mismatch = 0)
			if pam: 
				for p in pam:
					if ((p - 26) < 0) or ((p+8)>len(self.sequence_gRNA)): continue
					pos = self.start_gRNA + p
					seq = self.sequence_gRNA[p-26:p+8]
					if len(seq) != 34: continue
					count += 1
					gRNA_temp = Primer(f'gRNA', f'{self.id}-gRNA_{count}', seq, pos, source = 'gRNA', reverse = True)
					self.gRNAs.append(gRNA_temp)

		#if self.gRNAs: print (f"FOUND SOME {len(self.gRNAs)} gRNAs")

	def _predict_gRNAs(self, deepCpf1):
		if not self.gRNAs: 
			return

		if not deepCpf1.state:
			deepCpf1.state = deepCpf1.check_function()

		if not deepCpf1.state: 
			## PRINT WARNING 
			logger.info("deepCpf1 Is not working optimally!")
		
		temp_gRNA = []
		for g in self.gRNAs:
			sequence = [g.seq]
			score = deepCpf1.scoreSeqs(sequence) 
			if score[0][0] > self.threshold:
				g.score = round(score[0][0],5)
				print (g.id, g.seq, g.reverse, g.seq[4:31], reverseC(g.seq[4:31]), g.score)
				g.seq = g.seq[4:31]
				g.len = gRNA_LENGTH
				temp_gRNA.append(g)

		self.gRNAs = temp_gRNA

	##---------------------------------------------------------------------------------------------
	## Graph for LAMP Primer

	@staticmethod
	def create_relative_postional_plot(LAMP, prefix, conservation = [], entropy = []):
		os.chdir(prefix + PLOT_DIR)
		LAMP_title = LAMP.id

		xlim_start = LAMP.start 
		xlim_end = LAMP.end
		start_inner = LAMP.start_gRNA
		end_inner = LAMP.end_gRNA

		plt.switch_backend('agg')
		three = ['F3', 'B3']
		two = ['F2','B2']
		one = ['F1c', 'B1c']
		loop = ['LB', 'LF']

		if conservation and entropy: 
			#ax1 = plt.subplot(311) 
			fig, ax = plt.subplots(2,1,figsize=(15, 5), sharex=True)
		else: 
			fig, ax = plt.subplots(1,1,figsize=(15, 5), sharex=True)

		for primer in PRIMER_LIST:
			current_primer = getattr(LAMP, primer)
			if not current_primer: continue
			colour = 'k'
			primer_region = None
			if primer in three:
				primer_region = 4
			elif primer in two:
				primer_region = 3
			elif primer in one:
				primer_region = 1
			elif primer in loop:
				primer_region = 2
				colour = 'b'

			temp_dict = {}
			for i in range(current_primer.pos, current_primer.pos + current_primer.len):
				temp_dict[i] = primer_region
			ax[0].plot(list(temp_dict.keys()), list(temp_dict.values()), color=colour, linewidth=2)
		
		counter = 1
		if LAMP.gRNAs: 
			for guide_RNA in LAMP.gRNAs:
				counter += 0.25
				temp_gRNA_dict = {}
				if not guide_RNA.reverse:
					for j in range(guide_RNA.consensus_pos[0], guide_RNA.consensus_pos[1]):
						temp_gRNA_dict[j] = counter 
				else:
					for j in range(guide_RNA.consensus_pos[1], guide_RNA.consensus_pos[0]):
						temp_gRNA_dict[j] = counter 
				ax[0].plot(list(temp_gRNA_dict.keys()), list(temp_gRNA_dict.values()), color='r', linewidth=2)

		#ax[0].xlim(xlim_start, xlim_end-1)
		ax[0].axis(xmin = xlim_start, xmax = xlim_end)
		x_range = [i for i in range(xlim_start, xlim_end)]
		
		#plt.ylabel('Primer', fontsize=6)
		if LAMP.gRNAs:
			ax[0].axvspan(start_inner, end_inner, color='r', alpha=0.3)
			custom_lines = [Line2D([0], [0], color='r', lw=2),
							Line2D([0], [0], color='k', lw=2),
							Line2D([0], [0], color='b', lw=2)]
			ax[0].legend(custom_lines, ['gRNA', 'LAMP Primers', 'LOOP regions'])
			y_ticks = ['Primer 1', 'Loop', 'Primer 2', 'Primer 3']
			ax[0].set_yticks([1,2,3,4], y_ticks)

		elif not LAMP.gRNAs:
			custom_lines = [Line2D([0], [0], color='k', lw=2),
							Line2D([0], [0], color='b', lw=2)]
			ax[0].legend(custom_lines, ['LAMP Primers', 'LOOP regions'])
			y_ticks = ['Primer 1', 'Loop', 'Primer 2', 'Primer 3']
			ax[0].set_yticks([1,2,3,4], y_ticks)

		if conservation and entropy:
			ax[1].plot(x_range, conservation, linestyle="", color="blue", marker="o", markersize = 1.2)
			ax[1].set_ylabel("Conservation %",color="blue",fontsize=14)
			ax[1].tick_params(axis="x", labelsize=12) 
			ax[1].tick_params(axis='x', which='both', labelsize=10, labelbottom=True)
			
			ax2 = ax[1].twinx()
			ax2.plot(x_range, entropy, linestyle="", color="red", marker="o", markersize = 1.2)
			ax2.set_ylabel("Entropy Score",color="red",fontsize=14)
			fig.suptitle('Entropy | Conservation | Positional plot for ' + LAMP_title)
		else:
			fig.suptitle('Positional plot for ' + LAMP_title)

		ax[0].set_xlabel('Nucleotide position')
		figure = plt.gcf() # get current figure
		figure.set_size_inches(12, 6)
		plt.savefig(prefix + '_' + LAMP_title + '_visualisation.png', dpi = 128)
		os.chdir("..")

##-------------------------------------------------------------------------------------------------
## Primer or gRNA class

class Primer():

	def __init__(self, primer_type, primer_id, seq, pos, source = None, reverse = False):
		self.primer_type = primer_type
		self.id = primer_id
		self.seq = reverseC(seq.upper()) if reverse else seq.upper() 
		self.source = source ##GLAPD/LAVA/OTHERS/gRNA
		self.len = len(seq)
		self.reverse = reverse ##Relative to Reference Sequence ## LAMP_PRIMER_REVERSE
		self.pos = int(pos)

		## gRNA
		self.score = 0

		## For Consensus Sequence Search
		self.mode = None
		self.consensus_query = ''
		self._consensus_pos = (-1, -1)
		self.found = False
		self.cigar = ''
		self.Intolerant = False
		self.mismatch = float('inf') 
		self.alternative = []

		self.checked = False
		self.c_mean = 0
		self.c_SD = 0
		self.e_sum = 0
		self.c_score_list = []
		self.e_score_list = []

		self.c_score_bowtie = 0
		self.offtargets = []

	def scoring(self, c_score_list, e_score_list, store = False):
		self.c_mean = round(np.mean(c_score_list), 5)
		conservation_variance = np.sum(((c_score_list - self.c_mean)**2) / (self.len - 1))
		self.c_SD = round(np.sqrt(conservation_variance),5)

		self.e_sum = round(np.sum(e_score_list), 5)
		
		#self.e_mean = round(np.mean(e_score_list), 5)
		#entropy_variance = np.sum(((e_score_list - self.e_mean)**2) / (self.len - 1))
		#self.e_SD = round(np.sqrt(entropy_variance),5)

		if store:
			self.c_score_list = c_score_list
			self.e_score_list = e_score_list

		self.checked = True

	@property
	def consensus_pos(self):
		return self._consensus_pos

	@consensus_pos.setter
	def consensus_pos(self, start):
		end = start + self.len
		self._consensus_pos = (start, end) #Relative to Strand

##-------------------------------------------------------------------------------------------------
## Opening LAVA | GLAPD Files

def parse_LAVA_lists(lava_csv):
	import re
	primer_dict = dict()
	open_csv = open(lava_csv, 'r')
	lava_dict = None
	
	position_dict = defaultdict(list)
	with open(lava_csv, 'r') as open_csv:
		position_list = []
		Lamp_dict = dict()
		line = open_csv.readlines()
		lamp_number = 0
		for i in range(0, len(line)):
			if line[i].startswith('>'):
				sequence = line[i+1]
				sequence = sequence.strip('\n')
				line_ = re.findall(r'\w+', line[i]) 
				primer = line_[1]

				#lamp_number = int(line_[0])

				if primer in LAVA_LOOP:
					primer = primer.replace('OOP','')
				
				Lamp_dict[primer] = sequence
				if primer == 'BL':
					position_dict[lamp_number] = position_list
					primer_dict[lamp_number] = Lamp_dict
					lamp_number += 1
					position_list = []
					Lamp_dict = dict()

				if 'locations' in line_:
					location_index = line_.index('locations')
					for i in range(location_index+1, len(line_)-1,2):
						start = line_[i]
						end = line_[i+1]
						position_list.append((int(start), int(end)))
	
	# primer_dict => Lamp_dict[primer] = sequence
	# position_dict => position_list.append((int(start), int(end)))
	return primer_dict, position_dict

def check_LAVA_seq(position_dict, reference_sequence, primer_dict):
	final_primer_dict = dict()
	
	for index, (key, value) in enumerate(position_dict.items()):
		temp_primer_dict = defaultdict(list)
		sequence_list = primer_dict[index].values()
		sequence_list = [x.upper() for x in sequence_list]
		for i, k in enumerate(value):
			start = k[0]
			end = k[1]
			sequence = reference_sequence[start:end+1].upper()
			if 'N' in sequence: 
				logger.info("{}, {} not found for lamp primer number {}".format(PRIMER_LIST[i], sequence, str(index)))
				break
			reverse_sequence = reverseC(sequence)

			#check 
			sequence_checked = False
			if sequence in sequence_list:
				sequence_checked = True
			elif reverse_sequence in sequence_list:
				sequence_checked = True
				sequence = reverse_sequence
			else:
				for seq in sequence_list:
					if sequence in seq:
						sequence_checked = True 
						break
					elif reverse_sequence in seq:
						sequence_checked = True
						sequence = reverse_sequence
						break
			
			if not sequence_checked:
				logger.info("{}, {} not found for lamp primer number {}".format(PRIMER_LIST[i], sequence, str(index)))
				break
			
			temp_primer_dict[PRIMER_LIST[i]] = [int(start), sequence]
			if i == 7:
				final_primer_dict['LAMP_' + str(index + 1) + '_LAVA'] = temp_primer_dict
				temp_primer_dict = defaultdict(list)
	
	return final_primer_dict

def parse_GLAPD_lists(glapd_csv):
	primer_dict = dict()
	open_csv = open(glapd_csv, 'r')
	Lamp_dict = None
	counter = 0
	for line in open_csv:
		line = line.strip('\n')
		if line.startswith('The'):
			if Lamp_dict:
				primer_dict['LAMP_' + str(counter) + '_GLAPD'] = Lamp_dict
			Lamp_dict = defaultdict(list)
			counter += 1
		else: 
			line = line.replace(' ','')
			line_ = line.split(',')
			if line.split(':')[-1] == 'NULL':
				continue
			sequence = line.split(':')[-1]

			sequence = sequence.upper()
			front = line_[0].split(':')
			region = front[0]
			pos = int(front[2])
			Lamp_dict[region] = [pos, sequence]
			primer_dict['LAMP_' + str(counter) + '_GLAPD'] = Lamp_dict
	
	primer_dict['LAMP_' + str(counter) + '_GLAPD'] = Lamp_dict
	open_csv.close()
	return primer_dict

def convert_dict_2_class(primer_dict, overlap = False, threshold = 0):
	primer_class = []
	for key, value in primer_dict.items():
		temp_set = LAMP_set(key, overlap, threshold = threshold)
		source = key.split('_')[-1]
		for k, v in value.items():
			reverse = k in PRIMER_LIST_REVERSE
			primer = Primer(k, f'{key}-{k}', v[1], v[0], source = source, reverse = reverse)
			#print (k, f'{key}-{k}', v[1], v[0], source, reverse)
			setattr(temp_set, k, primer)
		primer_class.append(temp_set)

	if not primer_class:
		raise ExecutionError("No LAMP found!?")
	return primer_class

def open_LAMP_files(args_stored, ref, primer_class):
	overlap, threshold, LAMP_csv = args_stored
	ref = ref.ReconstructSequence()
	file_type = None
	primer_dict = {}
	for i in list(LAMP_csv): #process more than one file at once
		with open(i, 'r') as f:
			for line in f:
				line = line.strip('\n')
				if line.startswith('The'):
					file_type = 'GLAPD'
					break
				elif line.startswith('>'):
					file_type = 'LAVA'
					break
				else:
					raise ExecutionError('File type is wrong! Please enter output file from LAVA or GLAPD')
			
		if file_type == 'GLAPD':
			primer_dict_temp = parse_GLAPD_lists(i)
		elif file_type == 'LAVA':
			temp_primer_dict, position_dict = parse_LAVA_lists(i)
			primer_dict_temp = check_LAVA_seq(position_dict, ref, temp_primer_dict)
			
		if primer_dict:
			counter = len(primer_dict)
			primer_dict_temp = {'LAMP_' + str(int(k.split('_')[1]) + counter) + '_' + file_type:v for k,v in primer_dict_temp.items()}
			primer_dict.update(primer_dict_temp)
		else:
			primer_dict = primer_dict_temp


	primer_class += convert_dict_2_class(primer_dict, overlap, threshold)
	return primer_class

##-------------------------------------------------------------------------------------------------
## Multiprocessing

def partition(lst, n):
	division = len(lst) / n
	return [lst[round(division * i):round(division * (i + 1))] for i in range(n)]

def ranges(N, nb):
	temp_list = [(r.start, r.stop) for r in partition(range(len(N)), nb)]
	for i in temp_list: yield (i[0],i[1])

def _update_primer(g, matches, entropy_list, conservation_list, store = False):
	if matches:
		best = matches.pop(0)
		g.consensus_pos = best.pos[0]
		g.consensus_query = reverseC(best.query) if g.reverse else best.query
		g.mismatch = best.editDistance
		#g.Intolerant = best.intolerant
		g.cigar = best.cigar 

		start, end = g.consensus_pos[0], g.consensus_pos[1]
		e_score_list = entropy_list[start:end]
		c_score_list = conservation_list[start:end]
		g.scoring(c_score_list, e_score_list, store)
		g.found = True

		#print (1, g.id, g.c_mean, g.c_SD)
		if matches: g.alternative += matches
				
	else:
		g.found = False

def process_LAMP(args, 
				 consensus, 
				 reference, 
				 LAMP_list, 
				 queue, 
				 deepCpf1,
				 conservation_list, 
				 entropy_list):
	
	mismatch = args.mismatch
	store = True if args.keep_temp else False

	## Find gRNAs
	if args.find_gRNA:
		reference = reference.ReconstructSequence()
		for LAMP in LAMP_list:
			LAMP.get_Cas12gRNAs(reference, deepCpf1)
	del reference
	
	## Find MSA
	updated_LAMP = []
	if conservation_list and entropy_list:
		consensus_length = consensus.length
		consensus = consensus.ReconstructSequence()
		for LAMP in LAMP_list:
			for p in PRIMER_LIST:
				primer = getattr(LAMP, p)
				if not primer: continue
				query = primer.seq 
				matches = bitapSearch(consensus, query, reverse = primer.reverse, mismatch = mismatch)
				_update_primer(primer, matches, entropy_list, conservation_list, store)
				#print (2, primer.id, primer.c_mean, primer.c_SD, primer.found, primer.mismatch)

			for g in LAMP.gRNAs:
				query = g.seq if not g.reverse else reverseC(g.seq)
				matches = bitapSearch(consensus, query, reverse = g.reverse, mismatch = mismatch)
				_update_primer(g, matches, entropy_list, conservation_list, store)
				#print (3, g.id, g.c_mean, g.c_SD, g.found, g.mismatch, g.score)

			LAMP._update_start_end()
			if args.plot:
				start_temp = LAMP.start
				end_temp = LAMP.end  # + WINDOW_SIZE if LAMP.end+WINDOW_SIZE < consensus_length else LAMP.end

				if start_temp < 0 or end_temp < 0: continue
				#print (f'Slicing {start_temp} , {end_temp}')
				e_score_list = entropy_list[start_temp:end_temp]
				c_score_list = conservation_list[start_temp:end_temp]
				#if LAMP.end+WINDOW_SIZE < consensus_length:
				#	e_score_list = get_moving_average(e_score_list, window_size = WINDOW_SIZE)
				#	c_score_list = get_moving_average(c_score_list, window_size = WINDOW_SIZE)
				LAMP_set.create_relative_postional_plot(LAMP, args.prefix, c_score_list, e_score_list)

			updated_LAMP.append(LAMP)

	queue.put(updated_LAMP)

def queue_sam(q, updated_primers):
	while True:
		item = q.get()
		if item == 'Done':
			return
		else:
			for p in item:
				updated_primers.append(p)

def multi_search(main_args, wavelet, reference, LAMP_class, MSA_ed = False, scores = None):
	#if not MSA_ed or main_args.find_gRNA: 
	#	return LAMP_class

	if main_args.plot:
		if not os.path.isdir(main_args.prefix + PLOT_DIR): 
			os.mkdir(main_args.prefix + PLOT_DIR) 

	if main_args.find_gRNA:
		from LAMPgRNAtor.Cas12gRNAtor import gRNA, gRNA_Prediction

	from multiprocessing import Process, Queue, Manager, RawArray
	number_processors = main_args.CPU
	current_processes = []
	q = Queue()
	manager = Manager()
	updated_primers = manager.list()
	deepCpf1_object = gRNA_Prediction() if main_args.find_gRNA else None
 
	conservation_list, entropy_list = None, None
	
	if MSA_ed:
		import ctypes
		conservation_list, entropy_list = scores
		conservation_list = RawArray(ctypes.c_float, conservation_list)
		entropy_list = RawArray(ctypes.c_float, entropy_list)

	current_processes = []
	for i in ranges(LAMP_class, number_processors):
		current_LAMP = LAMP_class[i[0]:i[1]]
		pp = Process(target = process_LAMP, args=(main_args, wavelet, reference, current_LAMP, q, deepCpf1_object, conservation_list, entropy_list))
		current_processes.append(pp)

	for pp in current_processes:
		pp.start()

	sink = Process(target=queue_sam, args=(q, updated_primers))
	sink.start()

	for pp in current_processes:
		pp.join()

	q.put('Done')
	sink.join()

	return list(updated_primers)

def bowtie_to_summary(file, mode = "offtarget"):
	fasta_dict = {}
	file = open(file, '+r')
	seq_count = 0
	for line in file:
		line = line.strip('\n')
		if line.startswith('@SQ'):
			match = re.findall(r"\tSN:(\S+)", line)
			name = match[0] if len(match) else "unknown_reference"
			seq_count += 1
			continue
		if line.startswith('@'): continue
		cols = line.split("\t")
		sgrna_name = cols[0]
		mapping_status = int(cols[1])
		genome_name = cols[2]
		if mapping_status == 16 or mapping_status == 0:
			if mode == "conservation":
				if not sgrna_name in list(fasta_dict.keys()):
					fasta_dict[sgrna_name] = [genome_name]
				else:
					fasta_dict[sgrna_name].append(genome_name)
		
			if mode == "offtarget":
				if not sgrna_name in list(fasta_dict.keys()):
					fasta_dict[sgrna_name] = [genome_name]
				else:
					#print(fasta_dict[sgrna_name])
					fasta_dict[sgrna_name].append(genome_name)

	if mode == "conservation":
		fasta_dict = {k:round(len(list(set(v)))/seq_count * 100, 5) for k,v in fasta_dict.items()}

	return fasta_dict