#!/usr/bin/env python
import unittest
from nose.tools import assert_almost_equal, assert_equal

import os

from corrbin.score import read_contigs_file, parse_contig_description,\
    RocAxisFuns, ExperimentData
from probin import dna

class test_score(object):
    def setup(self):
        dna.DNA.generate_kmer_hash(3)

    def teardown(self):
        reload(dna)
    # testing function: read_contigs_file
    def test_read_contigs_file(self):
        cur_dir = os.path.dirname(__file__)
        file_name = os.path.join(cur_dir,"fixtures/generated_contigs_test.fna")
        open_file = open(file_name, 'r')
        contigs = read_contigs_file(open_file)
        assert_equal(len(contigs),29)
        assert_equal(contigs[0].id, 'Ehrlichia_canis_Jake_uid58071')
        assert_equal(contigs[0].contig_id, 1000)
        assert_equal(contigs[0].family, "Anaplasmataceae")


    def test_parse_contig_descritpion(self):
        id_string = "my_genome_name_uid123_100 my_family|my_genus|my_genus my_species"
        contig_id_hash = parse_contig_description(id_string)
        assert_equal("my_genome_name_uid123",contig_id_hash["genome"])
        assert_equal(100,contig_id_hash["contig_id"])
        assert_equal("my_family", contig_id_hash["family"])
        assert_equal("my_genus", contig_id_hash["genus"])
        assert_equal("my_genus my_species",contig_id_hash["species"])
                     
    def test_experiment_data(self):
        axis_funs = RocAxisFuns("included_contigs_ratio","precision")
        d = ExperimentData(axis_funs)
        cur_dir = os.path.dirname(__file__)
        input_file = os.path.join(cur_dir,"fixtures/score_mul_test.tsv")
        d.load_data_frame(input_file)
        d.standardize()
        d.classify("family")
        d.classify("genus")
        d.classify("genome")

        df = d.classification["family"]
        assert_equal(df[df.index == 1010].real_classif.values[0],False)
        assert_equal(list(df.real_classif.values).count(True),28)

        df = d.classification["genus"]
        assert_equal(df[df.index == 1010].real_classif.values[0],False)
        assert_equal(df[df.index == 1013].real_classif.values[0],False)
        assert_equal(df[df.index == 1027].real_classif.values[0],False)
        assert_equal(list(df.real_classif.values).count(True),26)


        df = d.classification["genome"]
        assert_equal(df[df.index == 1000].real_classif.values[0],False)
        assert_equal(df[df.index == 1001].real_classif.values[0],False)
        assert_equal(df[df.index == 1003].real_classif.values[0],False)
        assert_equal(df[df.index == 1004].real_classif.values[0],False)
        assert_equal(df[df.index == 1007].real_classif.values[0],False)
        assert_equal(df[df.index == 1010].real_classif.values[0],False)
        assert_equal(df[df.index == 1013].real_classif.values[0],False)
        assert_equal(df[df.index == 1018].real_classif.values[0],False)
        assert_equal(df[df.index == 1027].real_classif.values[0],False)

        assert_equal(list(df.real_classif.values).count(True),20)

        
