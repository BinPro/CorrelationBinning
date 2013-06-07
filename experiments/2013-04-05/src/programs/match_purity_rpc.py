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

def main(input_file, output_file,rpc_file):
    contig_purity_df = p.io.parsers.read_table(input_file,sep='\t', index_col=0)
    rpc_df = p.io.parsers.read_table(rpc_file,sep='\t',index_col=0,header=None)
    
#    contig_id_ix_dict = {}
#    for ix in contig_purity_df.index:
#        contig_id = contig_purity_df.ix[ix,:]['contig']
#        contig_id_ix_dict[contig_id]=ix
#
#    contig_reads_dict = {}
#    for contig_id in rpc_df.index:
#        ix = contig_id_ix_dict[contig_id]
#        contig_reads_dict[ix] = rpc_df.ix[contig_id,1]

    contig_purity_df['tot_nr_reads'] = p.Series(rpc_df[1].values,index=rpc_df.index)
    contig_purity_df.to_csv(output_file,sep='\t')


if __name__=="__main__":
    parser = ArgumentParser()
    parser.add_argument('file',
                        help='specify input file, default is stdin')
    parser.add_argument('-o', '--output', 
                        help='specify the base name of the output file.  The default is stdout')
    parser.add_argument('--rpc',
                        help='reads per contig file')
    args = parser.parse_args()

    main(args.file,args.output,args.rpc)
