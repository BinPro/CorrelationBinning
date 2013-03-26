#!/usr/bin/env python
import unittest
from nose.tools import assert_almost_equal, assert_equal,\
    assert_is_none
import os
file_path = os.path.realpath(__file__)
program_path = os.path.abspath(os.path.join(file_path,"..","..","src/programs/"))
import sys
sys.path.append(program_path)
import score_contigs_dirichlet
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

class TestScoreContigsDirichlet(object):
    def setUp(self):
        pass
    def tearDown(self):
        reload(dna)
        reload(score_contigs_dirichlet)
    # testing function: main
    def test_main(self):
        cur_dir = os.path.dirname(__file__)
        parsed_file_name = os.path.join(cur_dir,"fixtures/parsed_gen_2_2_test.txt")
        contig_file_name = os.path.join(cur_dir,"fixtures/contigs_2_2_100_10000_start_pos_test.fna")
        taxonomy_input_file = open(parsed_file_name,'r')
        contig_input_file = open(contig_file_name,'r')
        kmer_length = 3
        contig_length = 1000
        dir_path = os.path.join(cur_dir,"fixtures/reference_genomes")
        with tempfile.NamedTemporaryFile() as tmp_file:
            with RedirectStdStreams(stdout=tmp_file):
                score_contigs_dirichlet.main(contig_input_file,taxonomy_input_file, dir_path,kmer_length,contig_length)
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

    def test_main_avoid_overfit(self):
        cur_dir = os.path.dirname(__file__)
        parsed_file_name = os.path.join(cur_dir,"fixtures/parsed_gen_only_2.txt")
        contig_file_name = os.path.join(cur_dir,"fixtures/identical_contigs_different_genomes.fna")
        taxonomy_input_file = open(parsed_file_name,'r')
        contig_input_file = open(contig_file_name,'r')
        kmer_length = 3
        contig_length = 10000
        dir_path = os.path.join(cur_dir,"fixtures/reference_genomes")
        with tempfile.NamedTemporaryFile() as tmp_file:
            with RedirectStdStreams(stdout=tmp_file):
                score_contigs_dirichlet.main(contig_input_file,taxonomy_input_file, dir_path,kmer_length,contig_length)
            tmp_file.seek(0)
            num_lines = sum(1 for line in tmp_file)
            # #contigs * #genomes + #header
            # 2*2+1=5
            assert_equal(num_lines,5)
            tmp_file.seek(0)
            header_line = tmp_file.readline()
            all_entries = tmp_file.read()
            entries = all_entries.strip().split("\n")
            scores = [line.split("\t")[0] for line in entries]
            first_contig_first_genome = scores[0]
            first_contig_second_genome = scores[2]
            second_contig_first_genome = scores[1]
            second_contig_second_genome = scores[3]
            contig_id_differ1 = float(first_contig_first_genome) - float(second_contig_first_genome)
            assert_equal(abs(contig_id_differ1)>0.1,True)
            contig_id_differ2 = float(first_contig_second_genome) - float(second_contig_second_genome)
            assert_equal(abs(contig_id_differ2)>0.1,True)

        
