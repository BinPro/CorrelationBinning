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

def main(input_file, output_file,read_length):
    df = p.io.parsers.read_table(input_file,sep='\t')
    
    contig_reads_dict = {}
    for contig_id,group in df.groupby('contig_id'):
        group['per_base_read'] = (group['pos_stop'].astype(int)-group['pos_start'].astype(int))*group['cov_depth'].astype(int)
        contig_reads_dict[contig_id]=group['per_base_read'].sum()/float(read_length)
    
    output_series = p.Series(contig_reads_dict.values(),
                             index=contig_reads_dict.keys())
    output_series.to_csv(output_file,sep='\t')


if __name__=="__main__":
    parser = ArgumentParser()
    parser.add_argument('file',
                        help='specify input file, default is stdin')
    parser.add_argument('-o', '--output', 
                        help='specify the base name of the output file.  The default is stdout')
    parser.add_argument('--read_length',default=100,
                        help='speciy the read length used')
    args = parser.parse_args()

    main(args.file,args.output,args.read_length)
