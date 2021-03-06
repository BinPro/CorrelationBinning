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

def main(input_file, first_data, last_data,output_file, number_of_genomes):
    full_df = p.io.parsers.read_table(input_file,sep='\t')
    columns = full_df.ix[:,first_data:last_data].columns
    column_sums = np.sum(full_df.ix[:,first_data:last_data].values, axis=0)
    for index, column_name in enumerate(columns):
        full_df[column_name]/=float(column_sums[index])

    totals = p.Series(np.sum(full_df.ix[:,first_data:last_data].values,axis=1),index=full_df.index)
    totals.sort()
    included_indices = list(totals[-number_of_genomes:].index)
    
    shuffle(included_indices)
    filtered_df = full_df.ix[included_indices]
    filtered_df.to_csv(output_file,sep='\t')


if __name__=="__main__":
    parser = ArgumentParser()
    parser.add_argument('file',
                        help='specify input file, default is stdin')
    parser.add_argument('--first_data', 
                        help='specify the first data column')
    parser.add_argument('--last_data', 
                        help='specify the last data column')
    parser.add_argument('-o', '--output', 
                        help='specify the base name of the output file.  The default is stdout')
    parser.add_argument('--number_of_genomes', type=int,
                        help='specify the number of data rows to be extracted')
    args = parser.parse_args()

    main(args.file, args.first_data, args.last_data,args.output, args.number_of_genomes)
