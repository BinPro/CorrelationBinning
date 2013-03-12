#!/usr/bin/env python
import unittest
from nose.tools import assert_almost_equal, assert_equal,\
    assert_is_none
import os
file_path = os.path.realpath(__file__)
program_path = os.path.abspath(os.path.join(file_path,"..","..","src/programs/"))
import sys
sys.path.append(program_path)
import matrix_row
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

class TestMatrixRow(object):
    def setUp(self):
        pass
    def tearDown(self):
        pass
    # testing function: main
    def test_main(self):
        cur_dir = os.path.dirname(__file__)
        first_scores_file_name = os.path.join(cur_dir,"fixtures/score_mul_test.tsv")
        second_scores_file_name = os.path.join(cur_dir,"fixtures/score_mul_test2.tsv")
        with tempfile.NamedTemporaryFile() as tmp_file:
            with RedirectStdStreams(stdout=tmp_file):
                output = matrix_row.main(\
                    [first_scores_file_name,\
                         second_scores_file_name],
                    "genome",
                    "precision")
                assert_is_none(output)
                tmp_file.seek(0)
                first_line = tmp_file.readline()
                output_l = first_line.strip().split(" & ")
        assert_almost_equal(float(output_l[0]),0.68965517)
        assert_almost_equal(float(output_l[1]),0.89655172)
