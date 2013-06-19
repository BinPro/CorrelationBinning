#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  5 16:34:19 2013

@author: binni
"""
import numpy as np
from argparse import ArgumentParser
from Bio import SeqIO
from probin.dna import DNA
import pandas as pd
from collections import defaultdict

parser = ArgumentParser(description="Clustering of metagenomic contigs")

parser.add_argument('file', 
    help='specify input file on FASTA format')
parser.add_argument('-k', "--kmer", default=4, type=int,
    help="the kmer length")

args = parser.parse_args()

DNA.generate_kmer_hash(args.kmer)

seqs = list(SeqIO.parse(args.file,"fasta"))
contigs = [DNA(seq.id,seq.seq.tostring().upper(),calc_sign=True) for seq in seqs]

for contig in contigs:
    contigs_np = np.zeros((DNA.kmer_hash_count,),dtype=int)
    for i,value in contig.signature.iteritems():
        contigs_np[i] = value
    print "{0},{1}".format(contig.id,",".join(map(str,contigs_np)))

