#!/usr/bin/env python
import unittest
from nose.tools import assert_almost_equal, assert_equal

import os

from corrbin.score import read_contigs_file, parse_contig_id

class test_score(object):
    # testing function: read_contigs_file
    def test_false_positive_rate1(self):
        cur_dir = os.path.dirname(__file__)
        file_name = os.path.join(cur_dir,"fixtures/generated_contigs_test.fna")
        open_file = open(file_name, 'r')
        contigs = read_contigs_file(open_file)
        assert_equal(len(contigs),203)


    def test_parse_contig_id(self):
        id_string = "my_genome_name_uid123_100 my_family|my_genus|my_specie"
        contig_id_hash = parse_contig_id(id_string)
        assert_equal("my_genome_name_uid123",contig_id_hash("genome"))
        assert_equal(100,contig_id_hash("contig_id"))
        assert_equal("my_family", contig_id_hash("family"))
        assert_equal("my_genus", contig_id_hash("genus"))
        assert_equal("my_specie",contig_id_hash("specie"))
                     
    


