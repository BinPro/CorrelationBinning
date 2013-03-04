#!/usr/bin/env python
import unittest
from nose.tools import assert_almost_equal, assert_equal

from corrbin.contig_generation import SampleSetting

# testing class: SampleSetting
def test_sample_setting():
    s_set = SampleSetting("genomes",1000,\
                              1100,1200,\
                              True)
    assert_equal(s_set.prio,"genomes")
    assert_equal(s_set.no_contigs,1000)
    assert_equal(s_set.contig_min_length,1100)
    assert_equal(s_set.contig_max_length,1200)
    assert_equal(s_set.debug_mode,True)

