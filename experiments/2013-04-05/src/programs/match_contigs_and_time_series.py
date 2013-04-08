#!/usr/bin/env python
from argparse import ArgumentParser
import fileinput
import numpy as np
import pandas as p
from corrbin.score import ExperimentData, RocAxisFuns
from corrbin.classification import classify_bool
from corrbin.io import read_contigs_file, print_contigs_time_series
import matplotlib.pyplot as plt
import sys
import probin.dna as dna


def main(time_series_file, contig_file, output_file, contig_length, number_of_reads, first_data, last_data):
    dna.DNA.generate_kmer_hash(2)
    time_series_df = p.io.parsers.read_table(time_series_file,sep='\t')
    contigs = read_contigs_file(contig_file)

    # Match contigs with a real time_series
    ##  Count number of contigs per genome
    genomes = [c.id for c in contigs]
    genomes.sort()

    ##  Match time series index with
    ts_index_and_count = zip(time_series_df.index,[genomes.count(g) for g in genomes])
    multiplied_time_series = []
    for time_series_i,count in ts_index_and_count:
        multiplied_time_series += [time_series_i]*count

    # Simulate new time series for each contig, with the real as base.
    for contig,ts_i in zip(contigs,multiplied_time_series):
        ts = time_series_df.ix[ts_i]
        contig.otu = ts["# OTU"]
        contig.time_series = ts[first_data:last_data]
    
    
    sample_headers =time_series_df.ix[:,first_data:last_data].columns
    print_contigs_time_series(contigs,output_file,sample_headers)


if __name__=="__main__":
    parser = ArgumentParser()
    parser.add_argument('file',
                        help='specify input time series file, default is stdin')
    parser.add_argument('--contigs_file',
                        help='specify input contigs file')
    parser.add_argument('-o', '--output', 
                        help='specify the base name of the output file.  The default is stdout')
    parser.add_argument('--contig_length', type=int,
                        help='specify the length of contigs')
    parser.add_argument('--number_of_reads',type=int,
                        help='speciy total number of reads for a single sample')
    parser.add_argument('--first_data',
                        help='specify first data point for the samples')
    parser.add_argument('--last_data',
                        help='specify last data point for the samples')
    args = parser.parse_args()

    if args.output:
        output_file = open(args.output,'w+')
    else:
        output_file = sys.stdout
    main(args.file,args.contigs_file, output_file, args.contig_length, args.number_of_reads, args.first_data, args.last_data)

    output_file.close()
