#!/usr/bin/env python
import unittest
from nose.tools import assert_almost_equal, assert_equal
import os

from corrbin.contig_generation import SampleSetting,\
    sample_contig, read_parsed_taxonomy_file

from probin.dna import DNA
from Bio import SeqIO

DNA.generate_kmer_hash(3)
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

# testing function: sample_contig, in debug mode
def test_sample_contig():
    s_set = SampleSetting("genomes",1000,\
                              1100,1100,\
                              True)
    dir = os.path.dirname(__file__)
    dna_file = os.path.join(dir,"fixtures/8M_genome.fna")
    genome_seq = list(SeqIO.parse(dna_file, "fasta"))
    genome = DNA(id="test01", seq=str(genome_seq[0].seq))
    contig = sample_contig(genome,s_set,1)
    true_contig = str(genome_seq[0].seq[0:1100])
    assert_equal(contig.full_seq,true_contig)
    assert_equal(contig.id, 'test01 contig')

# testing function: sample_contig, without debug mode
def test_sample_contig2():
    s_set = SampleSetting("genomes",1000,\
                              1200,1200,\
                              False)
    dir = os.path.dirname(__file__)
    dna_file = os.path.join(dir,"fixtures/8M_genome.fna")
    genome_seq = list(SeqIO.parse(dna_file, "fasta"))
    genome = DNA(id="test01", seq=str(genome_seq[0].seq))
    contig = sample_contig(genome,s_set,1)
    assert_equal(len(contig.full_seq),1200)
    assert_equal(contig.id, 'test01 contig')
    
# testing function: read_parsed_taxonomy_file
def test_read_parsed_taxonomy_file():
    dir = os.path.dirname(__file__)
    file_name = os.path.join(dir,"fixtures/parsed_gen_2_2_test.txt")
    open_file = open(file_name, 'r')
    groups = read_parsed_taxonomy_file(open_file)
    assert_equal(len(groups),4)

