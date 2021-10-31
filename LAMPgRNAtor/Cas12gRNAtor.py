#!/usr/bin/env python3

# original source from https://github.com/MyungjaeSong/Paired-Library
# modified from Max Haeussler for better integration into parse LAMP 

from os.path import dirname, join
from numpy import *

import os
os.environ['KERAS_BACKEND'] = "theano"
os.environ['THEANO_FLAGS'] = "base_compiledir=temp/"
import logging
logging.getLogger('keras').setLevel(logging.WARNING)

from keras import backend as K
import importlib
importlib.reload(K)

from keras.models import Model
from keras.layers import Input, Conv1D
from keras.layers.merge import Multiply
from keras.layers.core import Dense, Dropout, Activation, Flatten
from keras.layers.convolutional import Convolution1D, AveragePooling1D

from LAMPgRNAtor.utils import PYTHON_FILE

def warn(*args, **kwargs):
	pass

import warnings
warnings.warn = warn


class gRNA():
	def __init__(self, start, end, prediction_score, forward = True):
		self.start = start
		self.end = end
		self.length = 27
		self.forward = forward
		self.prediction_score = None
		self.name = None

class gRNA_Prediction():

	def __init__(self):
		self.Seq_deepCpf1 = None
		self.DeepCpf1 = None

		self.create_model()
		self.state = self.check_function()

	def create_model(self):
		Seq_deepCpf1_Input_SEQ = Input(shape=(34,4))
		Seq_deepCpf1_C1 = Conv1D(80, 5, activation='relu')(Seq_deepCpf1_Input_SEQ)
		Seq_deepCpf1_P1 = AveragePooling1D(2)(Seq_deepCpf1_C1)
		Seq_deepCpf1_F = Flatten()(Seq_deepCpf1_P1)
		Seq_deepCpf1_DO1 = Dropout(0.3)(Seq_deepCpf1_F)
		Seq_deepCpf1_D1 = Dense(80, activation='relu')(Seq_deepCpf1_DO1)
		Seq_deepCpf1_DO2 = Dropout(0.3)(Seq_deepCpf1_D1)
		Seq_deepCpf1_D2 = Dense(40, activation='relu')(Seq_deepCpf1_DO2)
		Seq_deepCpf1_DO3 = Dropout(0.3)(Seq_deepCpf1_D2)
		Seq_deepCpf1_D3 = Dense(40, activation='relu')(Seq_deepCpf1_DO3)
		Seq_deepCpf1_DO4 = Dropout(0.3)(Seq_deepCpf1_D3)
		Seq_deepCpf1_Output = Dense(1, activation='linear')(Seq_deepCpf1_DO4)
		self.Seq_deepCpf1 = Model(inputs=[Seq_deepCpf1_Input_SEQ], outputs=[Seq_deepCpf1_Output])
		
		DeepCpf1_Input_SEQ = Input(shape=(34,4))
		DeepCpf1_C1 = Conv1D(80, 5, activation='relu')(DeepCpf1_Input_SEQ)
		DeepCpf1_P1 = AveragePooling1D(2)(DeepCpf1_C1)
		DeepCpf1_F = Flatten()(DeepCpf1_P1)
		DeepCpf1_DO1 = Dropout(0.3)(DeepCpf1_F)
		DeepCpf1_D1 = Dense(80, activation='relu')(DeepCpf1_DO1)
		DeepCpf1_DO2 = Dropout(0.3)(DeepCpf1_D1)
		DeepCpf1_D2 = Dense(40, activation='relu')(DeepCpf1_DO2)
		DeepCpf1_DO3 = Dropout(0.3)(DeepCpf1_D2)
		DeepCpf1_D3_SEQ = Dense(40, activation='relu')(DeepCpf1_DO3)
		
		DeepCpf1_Input_CA = Input(shape=(1,))
		DeepCpf1_D3_CA = Dense(40, activation='relu')(DeepCpf1_Input_CA)
		DeepCpf1_M = Multiply()([DeepCpf1_D3_SEQ, DeepCpf1_D3_CA])
		
		DeepCpf1_DO4 = Dropout(0.3)(DeepCpf1_M)
		DeepCpf1_Output = Dense(1, activation='linear')(DeepCpf1_DO4)
		self.DeepCpf1 = Model(inputs=[DeepCpf1_Input_SEQ, DeepCpf1_Input_CA], outputs=[DeepCpf1_Output])

		self.Seq_deepCpf1.load_weights(join(PYTHON_FILE, 'misc/Seq_deepCpf1_weights.h5'))
		self.DeepCpf1.load_weights(join(PYTHON_FILE, 'misc/DeepCpf1_weights.h5'))

	def scoreSeqs(self, seqs):
		for s in seqs:
			assert(len(s)==34)
		
		SEQ, CA0, CA1 = self.PREPROCESS(seqs)

		seqScores = self.Seq_deepCpf1.predict([SEQ], batch_size=50, verbose=0)
		deepScores0 = self.DeepCpf1.predict([SEQ, CA0], batch_size=50, verbose=0) * 3
		deepScores1 = self.DeepCpf1.predict([SEQ, CA1], batch_size=50, verbose=0) * 3
		
		#print (seqScores)
		seqScores = [x[0] for x in seqScores] 
		deepScores0 = [x[0] for x in deepScores0]
		deepScores1 = [x[0] for x in deepScores1]
		
		return seqScores, deepScores0, deepScores1

	def PREPROCESS(self, seqs):
		data_n = len(seqs)
		SEQ = zeros((data_n, 34, 4), dtype=int)
		CA0 = zeros((data_n, 1), dtype=int)
		CA1 = zeros((data_n, 1), dtype=int)
		
		for l in range(0, data_n):
			seq = seqs[l]
			for i in range(34):
				if seq[i] in "Aa":
					SEQ[l, i, 0] = 1
				elif seq[i] in "Cc":
					SEQ[l, i, 1] = 1
				elif seq[i] in "Gg":
					SEQ[l, i, 2] = 1
				elif seq[i] in "Tt":
					SEQ[l, i, 3] = 1

			CA0[l,0] = 0
			CA1[l,0] = 100

		return SEQ, CA0, CA1

	def check_function(self):
		sequence = ['TGACTTTGAATGGAGTCGTGAGCGCAAGAACGCT', 'GTTATTTGAGCAATGCCACTTAATAAACATGTAA'] 
		output = [55.699306, 53.469837]
		scoreSeqs_, _, _ = self.scoreSeqs(sequence)
		scoreSeqs_ = [round(float(i),6) for i in scoreSeqs_]
		return scoreSeqs_ == output
