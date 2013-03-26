#!/usr/bin/env python
import unittest
from nose.tools import assert_almost_equal, assert_equal,\
    assert_is_none
import os
import sys
file_path = os.path.realpath(__file__)
program_path = os.path.abspath(os.path.join(file_path,"..","..","src/programs/"))
sys.path.append(program_path)
import split_genomes
import tempfile

from probin import dna
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
        l = 10000
        with tempfile.NamedTemporaryFile() as tmp_file:
            with RedirectStdStreams(stdout=tmp_file):
                split_genomes.main(input_file,dir_path,l)
            tmp_file.seek(0)
            parts = list(SeqIO.parse(tmp_file, "fasta"))
            assert_equal(len(parts),1788)
            assert_equal(len(parts[0].seq),10000)
            genomes = [str(contig_seq.id).split("_")[0:-2]\
                           for contig_seq in parts]
            assert_equal(['Ehrlichia', 'canis', 'Jake', 'uid58071'],
                         genomes[0])
