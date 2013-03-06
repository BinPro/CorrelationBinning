#!/usr/bin/env python
import unittest
from nose.tools import assert_almost_equal, assert_equal,\
    assert_is_none
import os
import sys
sys.path.append('../src/programs/')
import generate_contigs
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

class TestGenerateContigs(object):
    def setUp(self):
        pass
    def tearDown(self):
        reload(dna)
    # testing function: main
    def test_main(self):
        cur_dir = os.path.dirname(__file__)
        parsed_file_name = os.path.join(cur_dir,"fixtures/parsed_gen_2_2_test.txt")
        dir_path = os.path.join(cur_dir,"fixtures/reference_genomes")
        input_file = open(parsed_file_name,'r')
        s_set = SampleSetting("genomes",10,\
                                  1100,1100,\
                                  True)
        with tempfile.NamedTemporaryFile() as tmp_file:
            with RedirectStdStreams(stdout=tmp_file):
                generate_contigs.main(input_file,dir_path,s_set)
            tmp_file.seek(0)
            contig_seqs = list(SeqIO.parse(tmp_file, "fasta"))
            assert_equal(len(contig_seqs),29)
            genomes = ["_".join(str(contig_seq.id).split("_")[0:-1])\
                           for contig_seq in contig_seqs]
            assert_equal("Capnocytophaga_canimorsus_Cc5_uid70727" \
                             in genomes,True)

