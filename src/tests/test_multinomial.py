#!/usr/bin/env python
import unittest
from nose.tools import assert_almost_equal, assert_equal,\
    assert_is_none
import os

from corrbin.multinomial import Experiment
from corrbin.contig_generation import SampleSetting

from probin import dna
from Bio import SeqIO

class TestMultinomial(object):
    def test_count_per_g(self):
        # testing class: Experiment, function count_per_g
        s_set = SampleSetting("genomes",10,\
                                  1100,1200,\
                                  True)
        test = Experiment(s_set, ["list","of","length","4"],[],None)
        assert_equal(test.count_per_g, [3,3,3,3])

    def test_count_per_g2(self):
        # testing class: Experiment, function count_per_g
        s_set = SampleSetting("genomes",10,\
                                  1100,1200,\
                                  True)
        test = Experiment(s_set, [5]*5, [], None)
        assert_equal(test.count_per_g,[2,2,2,2,2])

    def test_count_per_g3(self):
        # testing class: Experiment, function count_per_g
        s_set = SampleSetting("genomes",10,\
                                  1100,1200,\
                                  True)
        test = Experiment(s_set, [5]*5, [], None)
        assert_equal(test.count_per_g,[2,2,2,2,2])
        
    def test_count_per_g4(self):
        # testing class: Experiment, function count_per_g
        s_set = SampleSetting("groups",10,\
                                  1100,1200,\
                                  True)
        test = Experiment(s_set, [5]*4, [], None)
        assert_equal(test.count_per_g,[3,3,2,2])
        assert_equal(sum(test.count_per_g),10)
