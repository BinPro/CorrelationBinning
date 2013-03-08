#!/usr/bin/env python
import unittest
from nose.tools import assert_almost_equal, assert_equal,\
    assert_is_none
import os
file_path = os.path.realpath(__file__)
program_path = os.path.abspath(os.path.join(file_path,"..","..","src/programs/"))
import sys
sys.path.append(program_path)
import score_contigs_multinomial
import tempfile

from probin import dna
from corrbin.contig_generation import SampleSetting
from Bio import SeqIO

class RedirectStdStreams(object):
    def __init__(self, stdout=None, stderr=None):
        self._stdout = stdout or sys.stdout
        self._stderr = stderr or sys.stderr

    def __enter__(self):
        self.old_stdout, self.old_stderr = sys.stdout, sys.stderr
        self.old_stdout.flush(); self.old_stderr.flush()
        sys.stdout, sys.stderr = self._stdout, self._stderr

    def __exit__(self, exc_type, exc_value, traceback):
        self._stdout.flush(); self._stderr.flush()
        sys.stdout = self.old_stdout
        sys.stderr = self.old_stderr

class TestScoreContigsMultinomial(object):
    def setUp(self):
        pass
    def tearDown(self):
        reload(dna)
    # testing function: main
    def test_main(self):
        cur_dir = os.path.dirname(__file__)
        parsed_file_name = os.path.join(cur_dir,"fixtures/parsed_gen_2_2_test.txt")
        input_file = open(parsed_file_name,'r')
        kmer_length = 3
        dir_path = os.path.join(cur_dir,"fixtures/reference_genomes")
        s_set = SampleSetting("genomes",10,\
                                  1100,1100,\
                                  True)
        with tempfile.NamedTemporaryFile() as tmp_file:
            with RedirectStdStreams(stdout=tmp_file):
                score_contigs_multinomial.main(input_file,dir_path,kmer_length,s_set)
            tmp_file.seek(0)
            num_lines = sum(1 for line in tmp_file)
            # #contigs * #genomes + #header
            # 29*7+1=204
            assert_equal(num_lines,204)
            tmp_file.seek(0)
            header_line = tmp_file.readline()
            first_data_line = tmp_file.readline()
            elements =first_data_line.strip().split("\t")
            real_elements = ['Anaplasmataceae', 'Ehrlichia', 'Ehrlichia canis', 'Ehrlichia_canis_Jake_uid58071', 'Anaplasmataceae', 'Ehrlichia', 'Ehrlichia canis', 'Ehrlichia_canis_Jake_uid58071', '1000']

            assert_equal(real_elements, elements[1:])

