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


def main(time_series_file, contig_file, output_file, first_data, last_data,random_abundance=False, data_otu=False):
    dna.DNA.generate_kmer_hash(2)
    time_series_df = p.io.parsers.read_table(time_series_file,sep='\t')
    contig_df = p.io.parsers.read_table(contig_file, sep='\t', index_col = 0)

    # Match contigs with a real time_series
    ##  Count number of contigs per genome
    genomes = contig_df['full_read_mappings_strain'].unique()
    if "N/A" in genomes:
        sys.stderr.write("Please remove any 'N/A':s from the contig file" + '\n')
        sys.exit(-1)

    ## check that the number of genomes match the number of time series.
    if len(genomes) != len(time_series_df.index.values):
        sys.stderr.write("Number of genomes != number of time_series" + '\n')
        sys.exit(-1)
    
    ## Divide the number of reads for each genome with the total.
    total_number_of_reads_in_sample = np.sum(contig_df['full_read_mappings'])
    grouped_contig_df = contig_df.groupby('full_read_mappings_strain')

    otu_dict= {}
    for column in time_series_df.ix[:,first_data:last_data].columns:
        cr_dict = {}
        ## decide how many reads each genome gets
        if random_abundance:
            current_sample = np.random.multinomial(total_number_of_reads_in_sample,time_series_df[column].values)
        else:
            current_sample = time_series_df[column].values * total_number_of_reads_in_sample
        ## loop over this vector zipped with the contig group df
        for index_n,group in zip(enumerate(current_sample), grouped_contig_df):
            index = index_n[0]
            n = index_n[1]
            ## Spread the reads out over the contigs.
            group_sample = np.random.multinomial(int(n),group[1]['within_genome_read_ratio'].astype(float))
            ## save to some df.
            for contig_id,n_reads in zip(list(group[1].index),group_sample):
                cr_dict[contig_id] = n_reads
                if data_otu:
                    otu_dict[contig_id] = time_series_df['# OTU'].ix[index]
        contig_df[column] = p.Series(cr_dict.values(),index=cr_dict.keys())
    if data_otu:
        contig_df['data_otu'] = p.Series(otu_dict.values(),index=otu_dict.keys())
    contig_df.to_csv(output_file, sep='\t')
            

if __name__=="__main__":
    parser = ArgumentParser()
    parser.add_argument('file',
                        help='specify input time series file, default is stdin')
    parser.add_argument('--contigs_file',
                        help='specify input contigs file')
    parser.add_argument('-o', '--output', 
                        help='specify the base name of the output file.  The default is stdout')
    parser.add_argument('--first_data',
                        help='specify first data point for the samples')
    parser.add_argument('--last_data',
                        help='specify last data point for the samples')
    parser.add_argument('--random_abundance',action='store_true',
                        help='Add this tag if the species abundance should be random also')
    parser.add_argument('--data_otu',action='store_true',
                        help='Add this tag if the data_otu column is desired in output')
    args = parser.parse_args()

    if args.output:
        output_file = open(args.output,'w+')
    else:
        output_file = sys.stdout
    main(args.file,args.contigs_file, output_file,args.first_data, args.last_data,random_abundance=args.random_abundance, data_otu=args.data_otu)

    output_file.close()
