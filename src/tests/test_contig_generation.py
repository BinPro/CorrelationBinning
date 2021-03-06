#!/usr/bin/env python
import unittest
from nose.tools import assert_almost_equal, assert_equal,\
    assert_is_none
import os
import tempfile

from corrbin.contig_generation import SampleSetting,\
    sample_contig,SampleGroup 

from corrbin.misc import Uniq_id
from corrbin.io import read_parsed_taxonomy_file, read_FASTA_files
import probin.dna as dna
from Bio import SeqIO

class TestContigGeneration(object):
    def setUp(self):
        dna.DNA.generate_kmer_hash(3)
    def tearDown(self):
        reload(dna)
    # testing class: SampleSetting
    def test_sample_setting(self):
        s_set = SampleSetting("genomes",1000,\
                                  1100,1200,\
                                  True)
        assert_equal(s_set.prio,"genomes")
        assert_equal(s_set.no_contigs,1000)
        assert_equal(s_set.contig_min_length,1100)
        assert_equal(s_set.contig_max_length,1200)
        assert_equal(s_set.debug_mode,True)
        
    # testing function: sample_contig, in debug mode
    def test_sample_contig(self):
        s_set = SampleSetting("genomes",1000,\
                                  1100,1100,\
                                  True)
        dir = os.path.dirname(__file__)
        dna_file = os.path.join(dir,"fixtures/8M_genome.fna")
        genome_seq = list(SeqIO.parse(dna_file, "fasta"))
        genome = dna.DNA(id="test01", seq=str(genome_seq[0].seq))
        contig = sample_contig(genome,s_set,1)
        true_contig = str(genome_seq[0].seq[0:1100])
        assert_equal(contig.full_seq,true_contig)
        assert_equal(contig.id, 'test01 contig')
        
    # testing function: sample_contig, without debug mode
    def test_sample_contig2(self):
        s_set = SampleSetting("genomes",1000,\
                                  1200,1200,\
                                  False)
        dir = os.path.dirname(__file__)
        dna_file = os.path.join(dir,"fixtures/8M_genome.fna")
        genome_seq = list(SeqIO.parse(dna_file, "fasta"))
        genome = dna.DNA(id="test01", seq=str(genome_seq[0].seq))
        contig = sample_contig(genome,s_set,1)
        assert_equal(len(contig.full_seq),1200)
        assert_equal(contig.id, 'test01 contig')
        
    # Test function: generate_for_group
    def test_generate_group_contigs(self):
        cur_dir = os.path.dirname(__file__)
        parsed_file_name = os.path.join(cur_dir,"fixtures/parsed_gen_2_2_test.txt")
        open_file = open(parsed_file_name, 'r')
        groups = read_parsed_taxonomy_file(open_file)
        dir_path = os.path.join(cur_dir,"fixtures/reference_genomes")
        read_FASTA_files(groups, dir_path)
        uniq_id = Uniq_id(10)
        group = groups[-1]
        s_set = SampleSetting("genomes",20,\
                                  1100,1100,\
                                  True)
        sg = SampleGroup(s_set, group, uniq_id)
        sg.generate_group_contigs()
        
        assert_equal(len(sg.group.genomes[-1].contigs[-1].full_seq),1100)
        assert_equal(len(group.genomes[-1].contigs[-1].full_seq),1100)
        assert_equal(len(group.genomes[-1].contigs)+len(group.genomes[0].contigs), 20)


    def test_count_per_g(self):
        # testing class: SampleGroup, function count_per_g
        s_set = SampleSetting("genomes",10,\
                                  1100,1200,\
                                  True)
        sg = SampleGroup(s_set, ["list","of","length","4"],None)
        assert_equal(sg.count_per_g, [3,3,3,3])

    def test_count_per_g2(self):
        # testing class: SampleGroup, function count_per_g
        s_set = SampleSetting("genomes",10,\
                                  1100,1200,\
                                  True)
        sg = SampleGroup(s_set, [5]*5, None)
        assert_equal(sg.count_per_g,[2,2,2,2,2])

    def test_count_per_g3(self):
        # testing class: SampleGroup, function count_per_g
        s_set = SampleSetting("genomes",10,\
                                  1100,1200,\
                                  True)
        sg = SampleGroup(s_set, [5]*5, None)
        assert_equal(sg.count_per_g,[2,2,2,2,2])
        
    def test_count_per_g4(self):
        # testing class: SampleGroup, function count_per_g
        s_set = SampleSetting("groups",10,\
                                  1100,1200,\
                                  True)
        sg = SampleGroup(s_set, [5]*4, None)
        assert_equal(sg.count_per_g,[3,3,2,2])
        assert_equal(sum(sg.count_per_g),10)

