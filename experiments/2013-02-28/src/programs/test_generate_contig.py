# -*- coding: utf-8 -*-
#!/usr/bin/env python
from nose.tools import assert_almost_equal, assert_equal, assert_false, assert_in, assert_true
import generate_contigs
import tempfile
from os import linesep

class TestContigGeneration(object):
    TWO_GENOMES=""">genome1
GGGGCCCCTTTTTAAAATTATATGCGCGCGCAACACGG
>genome2
ATTATATATGAGAGCGCGCGCGGTGTGTCTCTGCTGC
"""
    OUTPUT_CONTIGS=""">genome1_0
GGGG
>genome1_1
CCCC
>genome1_2
TTTT
>genome1_3
TAAA
>genome1_4
ATTA
>genome1_5
TATG
>genome1_6
CGCG
>genome1_7
CGCA
>genome1_8
ACAC
>genome2_0
ATTA
>genome2_1
TATA
>genome2_2
TGAG
>genome2_3
AGCG
>genome2_4
CGCG
>genome2_5
CGGT
>genome2_6
GTGT
>genome2_7
CTCT
>genome2_8
GCTG
"""
       
    def test_get_two_genomes(self):
        with tempfile.NamedTemporaryFile() as tmp_file:
            tmp_file.write(self.TWO_GENOMES)
            tmp_file.seek(0)
            assert_equal(2,len(generate_contigs.get_sequences(tmp_file) ))
    
    def test_generate_contigs_of_length_4(self):
        with tempfile.NamedTemporaryFile() as tmp_file:
            tmp_file.write(self.TWO_GENOMES)
            tmp_file.seek(0)
            seqs = generate_contigs.get_sequences(tmp_file)
            contigs = generate_contigs.get_contigs(seqs,4)
            for k,v in contigs.iteritems():
                for i in range(len(v)):
                    v[i] = str(v[i])
            assert_equal({"genome1":["GGGG","CCCC","TTTT","TAAA","ATTA","TATG","CGCG","CGCA","ACAC"],  
                          "genome2":["ATTA","TATA","TGAG","AGCG","CGCG","CGGT","GTGT","CTCT","GCTG"]}, 
                          contigs)
    
    def test_write_to_fasta(self):
        with tempfile.NamedTemporaryFile() as tmp_file:
            tmp_file.write(self.TWO_GENOMES)
            tmp_file.seek(0)
            seqs = generate_contigs.get_sequences(tmp_file)
        contigs = generate_contigs.get_contigs(seqs,4)
        
        with tempfile.NamedTemporaryFile() as tmp_file:
            generate_contigs.write_sequences(contigs,tmp_file)
            tmp_file.seek(0)
            assert_equal(self.OUTPUT_CONTIGS,tmp_file.read())
            
            