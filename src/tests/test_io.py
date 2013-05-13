import unittest
from nose.tools import assert_almost_equal, assert_equal,\
    assert_is_none
import os 
import tempfile
from Bio import SeqIO

from corrbin.io import read_parsed_taxonomy_file, \
    genome_info_from_parsed_taxonomy_file, \
    read_FASTA_files, read_FASTA_files_no_groups, \
    print_parts

from corrbin.misc import Uniq_id
from corrbin.contig_generation import SampleSetting, \
    SampleGroup

from probin import dna

class TestIO(object):
    def setUp(self):
        dna.DNA.generate_kmer_hash(3)
    def tearDown(self):
        reload(dna)
    # testing function: read_parsed_taxonomy_file
    def test_read_parsed_taxonomy_file(self):
        cur_dir = os.path.dirname(__file__)
        file_name = os.path.join(cur_dir,"fixtures/parsed_gen_2_2_test.txt")
        open_file = open(file_name, 'r')
        groups = read_parsed_taxonomy_file(open_file)
        assert_equal(len(groups),3)
        
    # testing function: genomes_from_parsed_taxonomy_file
    def test_genome_info_from_parsed_taxonomy_file(self):
        cur_dir = os.path.dirname(__file__)
        file_name = os.path.join(cur_dir,"fixtures/parsed_gen_2_2_test.txt")
        open_file = open(file_name, 'r')
        genomes = genome_info_from_parsed_taxonomy_file(open_file)
        assert_equal(len(genomes),7)
        assert_equal(genomes[0]["genus"], "Ehrlichia")
        

    # Test function: read_FASTA_files
    def test_read_FASTA_files(self):
        cur_dir = os.path.dirname(__file__)
        parsed_file_name = os.path.join(cur_dir,"fixtures/parsed_gen_2_2_test.txt")
        open_file = open(parsed_file_name, 'r')
        groups = read_parsed_taxonomy_file(open_file)
        dir_path = os.path.join(cur_dir,"fixtures/reference_genomes")
        output = read_FASTA_files(groups, dir_path)
        assert_is_none(output)
        last_genome = groups[-1].genomes[-1]
        assert_equal(len(last_genome.full_seq),2612925)
        assert_equal(last_genome.id, "Capnocytophaga_ochracea_DSM_7271_uid59197")
        # Same family and genera within group
        assert_equal(groups[-1].genomes[-1].family, groups[-1].genomes[0].family)
        assert_equal(groups[-1].genomes[-1].genus, groups[-1].genomes[0].genus)
        # A correct family
        assert_equal(groups[-1].genomes[-1].family, "Flavobacteriaceae")


    def test_read_single_FASTA_file(self):
        cur_dir = os.path.dirname(__file__)
        parsed_file_name = os.path.join(cur_dir,"fixtures/parsed_gen_0_0_mock_test_complete.txt")
        open_file = open(parsed_file_name, 'r')
        groups = read_parsed_taxonomy_file(open_file)
        dir_path = os.path.join(cur_dir,"fixtures/mock_references_test.fa")
        output = read_FASTA_files(groups, dir_path, dir_structure='single_fasta_file')
        assert_is_none(output)
        last_genome = groups[-1].genomes[-1]
        assert_equal(len(last_genome.full_seq),2222430) #wrong length
        assert_equal(last_genome.id, "Pyrobaculum_aerophilum_str._IM2")
        # Same family and genera within group
        assert_equal(groups[-1].genomes[-1].family, groups[-1].genomes[0].family)
        assert_equal(groups[-1].genomes[-1].genus, groups[-1].genomes[0].genus)
        # A correct family
        assert_equal(groups[-1].genomes[-1].family, "Thermoproteaceae") # wrong family

    def test_read_signle_FASTA_file_no_groups(self):
        cur_dir = os.path.dirname(__file__)
        parsed_file_name = os.path.join(cur_dir,"fixtures/parsed_gen_0_0_mock_test_complete.txt")
        open_file = open(parsed_file_name, 'r')
        meta_genomes = genome_info_from_parsed_taxonomy_file(open_file)
        dir_path = os.path.join(cur_dir,"fixtures/mock_references_test.fa")
        real_genomes = read_FASTA_files_no_groups(meta_genomes, dir_path, dir_structure='single_fasta_file')
        assert_equal(len(real_genomes),15)
        last_genome = real_genomes[-1]
        assert_equal(len(last_genome.full_seq),2222430)
        assert_equal(last_genome.id, "Pyrobaculum_aerophilum_str._IM2")
        # A correct family
        assert_equal(real_genomes[-1].family, "Thermoproteaceae")

    def test_read_FASTA_files_no_groups(self):
        cur_dir = os.path.dirname(__file__)
        parsed_file_name = os.path.join(cur_dir,"fixtures/parsed_gen_2_2_test.txt")
        open_file = open(parsed_file_name, 'r')
        meta_genomes = genome_info_from_parsed_taxonomy_file(open_file)
        dir_path = os.path.join(cur_dir,"fixtures/reference_genomes")
        real_genomes = read_FASTA_files_no_groups(meta_genomes, dir_path)
        assert_equal(len(real_genomes),7)
        last_genome = real_genomes[-1]
        assert_equal(len(last_genome.full_seq),2612925)
        assert_equal(last_genome.id, "Capnocytophaga_ochracea_DSM_7271_uid59197")
        # A correct family
        assert_equal(real_genomes[-1].family, "Flavobacteriaceae")
    
    def test_print_group_contigs(self):
        cur_dir = os.path.dirname(__file__)
        parsed_file_name = os.path.join(cur_dir,"fixtures/parsed_gen_2_2_test.txt")
        open_file = open(parsed_file_name, 'r')
        groups = read_parsed_taxonomy_file(open_file)
        dir_path = os.path.join(cur_dir,"fixtures/reference_genomes")
        read_FASTA_files(groups, dir_path)
        uniq_id = Uniq_id(10)
        group = groups[-1]
        s_set = SampleSetting("genomes",10,\
                                  1100,1100,\
                                  True)
        sg = SampleGroup(s_set, group, uniq_id)
        sg.generate_group_contigs()
        
        with tempfile.NamedTemporaryFile() as tmp_file:
            sg.print_group_contigs(tmp_file)
            tmp_file.seek(0)
            contig_seqs = list(SeqIO.parse(tmp_file, "fasta"))
            assert_equal(len(contig_seqs),10)
            assert_equal(contig_seqs[0].id, "Capnocytophaga_canimorsus_Cc5_uid70727_10")
            d_string = "Capnocytophaga_canimorsus_Cc5_uid70727_10 Flavobacteriaceae|Capnocytophaga|Capnocytophaga canimorsus"
            assert_equal(contig_seqs[0].description, d_string)


    def test_print_parts(self):
        cur_dir = os.path.dirname(__file__)
        parsed_file_name = os.path.join(cur_dir,"fixtures/parsed_gen_2_2_test.txt")
        open_file = open(parsed_file_name, 'r')
        groups = read_parsed_taxonomy_file(open_file)
        dir_path = os.path.join(cur_dir,"fixtures/reference_genomes")
        read_FASTA_files(groups, dir_path)
        uniq_id = Uniq_id(10)
        with tempfile.NamedTemporaryFile() as tmp_file:
            for group_index, group in enumerate(groups):
                for genome in group.genomes:
                    parts = genome.split_seq(10000)
                    print_parts(parts,tmp_file,uniq_id, genome)
            tmp_file.seek(0)
            genome_parts = list(SeqIO.parse(tmp_file,"fasta"))
            assert_equal(len(genome_parts),1788)
            assert_equal(genome_parts[0].id,"Ehrlichia_canis_Jake_uid58071_10_0")
