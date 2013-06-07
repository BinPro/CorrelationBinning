#!/usr/bin/env python
from argparse import ArgumentParser
import fileinput
import numpy as np
from random import shuffle
import pandas as p
from corrbin.score import ExperimentData, RocAxisFuns
from corrbin.classification import classify_bool
import matplotlib.pyplot as plt
import sys

def main(input_file, read_mapping_column, contig_length_column, otu_column, output_file, read_length, sorted_output=False, remove_na=False):
    full_df = p.io.parsers.read_table(input_file,sep='\t',index_col=0)
    grouped = full_df.groupby(otu_column)
    cr_hash = {}
    read_ratio_hash = {}
    for otu,group in grouped:
        if otu != "N/A":
            read_mappings = eval("group." + read_mapping_column)
            total_read_mappings = np.sum(read_mappings.astype(float))
            for contig in group.index.values:
                read_ratio_hash[contig] = float(group.ix[contig][read_mapping_column])/total_read_mappings
                cr_hash[contig] = total_read_mappings
    read_series = p.Series(cr_hash.values(), index = cr_hash.keys())
    read_ratio_series = p.Series(read_ratio_hash.values(), index = read_ratio_hash.keys())

    full_df['tot_number_of_reads'] = read_series
    full_df['within_genome_read_ratio'] = read_ratio_series
    full_df['full_read_mappings'] = full_df[read_mapping_column]
    if sorted_output:
        full_df= full_df.sort_index(by=['tot_number_of_reads',otu_column])
        
    if remove_na:
        full_df = full_df[full_df[otu_column] != "N/A"]
        full_df = full_df.dropna()
    full_df.to_csv(output_file,sep='\t')


if __name__=="__main__":
    parser = ArgumentParser()
    parser.add_argument('file',
                        help='specify input file, default is stdin')
    parser.add_argument('--read_mapping_column',
                        help='specify the read mapping column')
    parser.add_argument('--contig_length_column',
                        help='specify the contig length column')
    parser.add_argument('--otu_column',
                        help='specify the otu column')
    parser.add_argument('--read_length', type=int,
                        help='specify the read length')
    parser.add_argument('-o', '--output', 
                        help='specify the base name of the output file.  The default is stdout')
    parser.add_argument('--sorted', action='store_true',
                        help='Add tag if output should be sorted by otu and otu average coverage')
    parser.add_argument('--remove_na', action='store_true',
                        help='Add tag if any "N/A":s should be removed from output')

    args = parser.parse_args()

    main(args.file, args.read_mapping_column, args.contig_length_column,args.otu_column, args.output, args.read_length, sorted_output= args.sorted, remove_na = args.remove_na)
