#!/usr/bin/env python
import unittest
from nose.tools import assert_almost_equal, assert_equal,\
    assert_is_none
import os
file_path = os.path.realpath(__file__)
program_path = os.path.abspath(os.path.join(file_path,"..","..","src/programs/"))
import sys
sys.path.append(program_path)
import score_time_series
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

class TestScoreTimeSeries(object):
    def setUp(self):
        pass
    def tearDown(self):
        reload(dna)
        reload(score_time_series)
    # testing function: main
    def test_main(self):
        cur_dir = os.path.dirname(__file__)
        file_names = ["fixtures/contigs_2_2_100_10000_start_pos_test.fna",
                      "fixtures/contigs_2_2_100_10000_start_pos_time_series_test.tsv",
                      "fixtures/time_series_test.tsv",
                      "fixtures/parsed_gen_2_2_test.txt",
                      "fixtures/reference_genomes"]
        files = tuple([open(os.path.join(cur_dir,file_name),'r')
                       for file_name in file_names[:-1]] + 
                      [os.path.join(cur_dir,file_names[-1])])
        contig_length = 10000
        total_read_count = 75000000
        assembly_length = 10**6
        with tempfile.NamedTemporaryFile() as tmp_file:
            with RedirectStdStreams(stdout=tmp_file):
                score_time_series.main(*files + 
                                        (contig_length,
                                         total_read_count, 
                                         assembly_length))
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

