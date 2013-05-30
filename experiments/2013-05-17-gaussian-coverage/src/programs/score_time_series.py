#!/usr/bin/env python
import fileinput
import sys
import os
import numpy as np
import pandas as p
from copy import copy
from argparse import ArgumentParser
from probin.model.coverage import isotropic_gaussian as model
from probin.model.coverage import log_coverage
from probin.dna import DNA
from Bio import SeqIO
from corrbin.misc import all_but_index, Uniq_id, GenomeGroup
from corrbin.multinomial import Experiment
from corrbin.score import Score
from corrbin.io import read_contigs_file, genome_info_from_parsed_taxonomy_file, read_FASTA_files_no_groups, read_time_series, read_time_series_file_genomes

def main(contigs_file,contig_time_series_file, genome_time_series_file, taxonomy_file,dir_path, kmer_length,first_data,last_data):

    DNA.generate_kmer_hash(kmer_length)

    contig_time_series_df = p.io.parsers.read_table(contig_time_series_file,sep='\t',index_col=0)

    include_contigs = {}
    for ix in contig_time_series_df.index:
       include_contigs[ix] = True

    contigs = read_contigs_file(contigs_file,start_position=False,taxonomy_info=False,filter_dict=include_contigs)
            
    contig_strain_otu_dic = {}

    strains = contig_time_series_df['full_read_mappings_strain'].values
    otus = contig_time_series_df['data_otu'].values
    for strain,otu in zip(strains,otus):
        contig_strain_otu_dic[strain] = otu

    total_number_of_reads_in_sample = np.sum(contig_time_series_df['full_read_mappings'])

    if len(contigs)!=len(contig_time_series_df.index):
        sys.stderr.write('Number of contigs: ' + str(len(contigs)) +'\n')
        sys.stderr.write('Number of time_series_indices: ' + str(len(contig_time_series_df.index))+'\n')
        
        raise TypeError("The number of contigs and time series does not match")


    for contig in contigs:
        contig.mapping_reads = contig_time_series_df.ix[contig.contig_id,first_data:last_data]
        contig.log_coverage = log_coverage.read_mappings_to_log_coverage(contig.mapping_reads,len(contig.full_seq),100)
        contig.species = contig_time_series_df.ix[contig.contig_id]['Contig_Species']
        contig.genus = contig_time_series_df.ix[contig.contig_id]['Contig_Genus']
        contig.family = contig_time_series_df.ix[contig.contig_id]['Contig_Family']

    # Divide genomes into groups, one for each genus
    meta_genomes = genome_info_from_parsed_taxonomy_file(taxonomy_file)

    # Fetch sequence for each genome
    genomes = read_FASTA_files_no_groups(meta_genomes, dir_path, dir_structure='single_fasta_file')

    # Fetch time series for each genome
    read_time_series_file_genomes(genomes, genome_time_series_file, contig_strain_otu_dic,first_data,last_data,total_number_of_reads_in_sample*100)

    for genome in genomes:
        mr = genome.mapping_reads
        if mr.shape[0] == 0:
            sys.stderr.write("skipping: " + str(genome.id)+'\n')
            continue
        else:
            sys.stderr.write("Including genome: " + str(genome.id) + '\n')
            sys.stderr.write("Time series: " + str(mr) + '\n')
        genome.pseudo_par = model.fit_nonzero_parameters(mr)

    scores = []
    for contig in contigs:
        for genome in genomes:
            p_val = model.log_pdf(\
                contig.log_coverage, *genome.pseudo_par)
            scores.append(\
                Score(p_val, contig, genome, contig.contig_id))

    sys.stdout.write("p_value\tcontig_family\tcontig_genus\tcontig_species\tcontig_genome\tcompare_family\tcompare_genus\tcompare_species\tcompare_genome\tcontig_id" + os.linesep)
    for score in scores:
        sys.stdout.write(str(score) + '\n')

if __name__=="__main__":
    parser = ArgumentParser()
    parser.add_argument('contigs', 
                        help='specify contig file with start positions, default is stdin')
    parser.add_argument('contig_time_series_file', 
                        help='Specify contig time series file')
    parser.add_argument('genome_time_series_file', 
                        help='Specify genome time series file')
    parser.add_argument('-o', '--output', 
                        help='specify the output file.  The default is stdout')
    parser.add_argument('-t', '--taxonomy',
                        help='specify the taxonomy file.')
    parser.add_argument('-k', '--kmer_length', default=4, type=int,
                        help='specify the kmer length, default is 4')
    parser.add_argument('-d', '--directory_path', default='/home/johannes/repos/DATA/reference_genomes_ncbi', 
                        type=str, help='specify the path to where the reference genomes are located locally')
    parser.add_argument('--first_data',
                        help='specify the name of the column containing the first abundance data')
    parser.add_argument('--last_data',
                        help='specify the name of the column containing the last abundance data')
    args = parser.parse_args()
    if args.output and args.output != '-':
        sys.stdout = open(args.output, 'w')
    if args.taxonomy:
        taxonomy_file = open(args.taxonomy, 'r')

    main(args.contigs,args.contig_time_series_file,args.genome_time_series_file,taxonomy_file, args.directory_path, args.kmer_length, args.first_data,args.last_data)

