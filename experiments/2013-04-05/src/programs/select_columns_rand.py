#!/usr/bin/env python
from argparse import ArgumentParser
import fileinput
import numpy as np
import pandas as p
from corrbin.score import ExperimentData, RocAxisFuns
from corrbin.classification import classify_bool
from corrbin.io import read_contigs_file, print_contigs_time_series
from corrbin.time_series import poisson as time_series_model
import matplotlib.pyplot as plt
import sys
import probin.dna as dna
import random

def main(time_series_file,output_file, first_data, last_data,nr_columns):
    ts_df = p.io.parsers.read_table(time_series_file,sep='\t',index_col=0)

    for i in nr_columns:
        while len(ts_df.ix[:,first_data:last_data].columns)>i:
            col = random.choice(ts_df.ix[:,first_data:last_data].columns)
            if col == first_data or col == last_data:
                continue
            else:
                del ts_df[col]
        ts_df.to_csv(output_file +"_" + str(i) + ".tsv", sep='\t')

            

if __name__=="__main__":
    parser = ArgumentParser()
    parser.add_argument('file',
                        help='specify input time series file, default is stdin')
    parser.add_argument('-o', '--output', 
                        help='specify the base name of the output file.  The default is stdout')
    parser.add_argument('--first_data',
                        help='specify first data point for the samples')
    parser.add_argument('--last_data',
                        help='specify last data point for the samples')
    parser.add_argument('--nr_columns', nargs='*',type=int, default=[46,30,15,2],
                        help='specify list of integers used as number of columns in the subsetted output files')
    args = parser.parse_args()

    main(args.file, args.output,args.first_data, args.last_data,args.nr_columns)
