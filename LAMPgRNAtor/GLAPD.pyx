#!/usr/bin/env python3

import os
from os.path import join as pjoin

PYTHON_FILE = os.path.dirname(os.path.abspath(__file__)) + '/LAMPgRNAtor'
MODULE_PATH = pjoin(PYTHON_FILE, 'src')
PAR_PATH = pjoin(MODULE_PATH, 'Primer3_Par') + '/'

cdef extern from "single.c":
	void SINGLE(const char* reference, const char* prefix, char* par_path, char* store_path)

cdef extern from "design.c":
	void DESIGNER(const char* input_py, const char* reference, int minLen_py, const char* prefix, char* par_path, char* store_path, int expect)

def single_candidate(reference: bytes, prefix: bytes, store_path: bytes) -> None:
	SINGLE(reference, prefix, str.encode(PAR_PATH), store_path)

def design_LAMP(prefix: bytes, reference: bytes, output: bytes, store_path: bytes, minlen: int = 0, expect: int = 100) -> None:
	DESIGNER(prefix, reference, minlen, output, str.encode(PAR_PATH), store_path, expect)

