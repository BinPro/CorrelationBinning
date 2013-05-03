#!/usr/bin/env python
import fileinput
import sys
import os
import numpy as np
from copy import copy
from argparse import ArgumentParser
from probin.model.coverage import binomial as model
from probin.dna import DNA
from Bio import SeqIO
from corrbin.misc import all_but_index, Uniq_id, GenomeGroup
from corrbin.multinomial import Experiment
from corrbin.score import Score
from corrbin.io import read_contigs_file, genome_info_from_parsed_taxonomy_file, read_FASTA_files_no_groups, read_time_series, read_time_series_file_genomes

def main(read_mapping_file_name,read_length):
    df_rm = p.io.parsers.read_table(read_mapping_file_name,sep='\t',index_col=0)
    cl = df_rm[df_rm.full_read_mappings != 'N/A'].contig_length
    frm= df_rm[df_rm.full_read_mappings != 'N/A'].full_read_mappings.astype(float)
    df_rm['coverage'] = frm*read_length/cl


if __name__=="__main__":
    parser = ArgumentParser()
    parser.add_argument('read_mapping_file', 
                        help='specify read mapping file with start positions, default is stdin')
    parser.add_argument('-o', '--output', 
                        help='specify the output file.  The default is stdout')
    parser.add_argument('-l', '--read_length', default=100, type=int,
                        help='specify the read length, default is 100')
    args = parser.parse_args()
    if args.output and args.output != '-':
        sys.stdout = open(args.output, 'w')
        
    main(args.read_mapping_file, args.read_length)

