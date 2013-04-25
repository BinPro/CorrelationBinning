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

def main(contigs_file,contig_time_series_file, genome_time_series_file, taxonomy_file,dir_path, contig_length, total_read_count,assembly_length,first_data,last_data):

    DNA.generate_kmer_hash(2)

    contigs = read_contigs_file(contigs_file,start_position=True)
    
    contig_time_series_df = read_time_series(contig_time_series_file)

    if len(contigs)!=len(contig_time_series_df.index):
        raise TypeError("The number of contigs and time series does not match")
    
    for contig in contigs:
        contig.mapping_reads = contig_time_series_df[contig_time_series_df.contig_id == contig.contig_id]

    # Divide genomes into groups, one for each genus
    meta_genomes = genome_info_from_parsed_taxonomy_file(taxonomy_file)

    # Fetch sequence for each genome
    genomes = read_FASTA_files_no_groups(meta_genomes, dir_path)

    # Fetch time series for each genome
    read_time_series_file_genomes(genomes, genome_time_series_file)

    for genome in genomes:
        genome.pseudo_par = model.fit_nonzero_parameters([genome],total_read_count)

    scores = []
    for contig in contigs:
        for genome in genomes:
            p_val = model.log_probability(\
                    contig, genome.pseudo_par, total_read_count,assembly_length)
            scores.append(\
                Score(p_val, contig, genome, contig.contig_id))

    sys.stdout.write("p_value\tcontig_family\tcontig_genus\tcontig_species\tcontig_genome\tcompare_family\tcompare_genus\tcompare_species\tcompare_genome\tcontig_id" + os.linesep)
    for score in scores:
        sys.stdout.write(str(score) + '\n')
   
if __name__=="__main__":
    parser = ArgumentParser()
    parser.add_argument('contigs', 
                        help='specify contig file with start positions, default is stdin')
    parser.add_argument('-o', '--output', 
                        help='specify the output file.  The default is stdout')
    parser.add_argument('-t', '--taxonomy',
                        help='specify the taxonomy file.')
    parser.add_argument('-k', '--kmer_length', default=4, type=int,
                        help='specify the kmer length, default is 4')
    parser.add_argument('-d', '--directory_path', default='/home/johannes/repos/DATA/reference_genomes_ncbi', 
                        type=str, help='specify the path to where the reference genomes are located locally')
    parser.add_argument('-c', '--contig_length', type=int,
                        help='specify the contig length to correctly avoid overfit')
    args = parser.parse_args()
    if args.output and args.output != '-':
        sys.stdout = open(args.output, 'w')
    if args.taxonomy:
        taxonomy_file = open(args.taxonomy, 'r')
        
    main(args.contigs, taxonomy_file, args.directory_path, args.kmer_length, args.contig_length)

