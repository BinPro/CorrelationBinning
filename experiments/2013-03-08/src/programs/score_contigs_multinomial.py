#!/usr/bin/env python
import fileinput
import sys
import os
from copy import copy,deepcopy
from argparse import ArgumentParser
from probin.model.composition import multinomial as mn 
from probin.dna import DNA
from Bio import SeqIO
from corrbin.misc import all_but_index, Uniq_id, GenomeGroup
from corrbin.multinomial import Experiment
from corrbin.contig_generation import SampleSetting 
from corrbin.score import Score
from corrbin.io import read_contigs_file, genome_info_from_parsed_taxonomy_file, read_FASTA_files_no_groups
import pandas as p

def main(contigs_file,taxonomy_file, dir_path, kmer_length,dir_structure, taxonomy_info_in_contigs, filter_file=None):

    groups = []
    DNA.generate_kmer_hash(kmer_length)
    if filter_file:
        filter_df = p.io.parsers.read_table(filter_file,sep='\t',index_col=0)
        include_contigs = filter_df.index
        filter_contig_dic = {}
        for contig in include_contigs:
            filter_contig_dic[contig] = True
            
        contigs = read_contigs_file(contigs_file,taxonomy_info=taxonomy_info_in_contigs,filter_dict=filter_contig_dic)
    else:
        contigs = read_contigs_file(contigs_file,taxonomy_info=taxonomy_info_in_contigs)

    # Divide genomes into groups, one for each genus
    meta_genomes = genome_info_from_parsed_taxonomy_file(taxonomy_file)

    # Fetch sequence for each genome
    genomes = read_FASTA_files_no_groups(meta_genomes, dir_path, dir_structure=dir_structure)

    for genome in genomes:
        genome.calculate_signature()
        genome.pseudo_par = mn.fit_nonzero_parameters([genome])

    scores = []
    for contig in contigs:
        if (not taxonomy_info_in_contigs) and (filter_file):
            contig.species = filter_df.ix[contig.contig_id]['Contig_Species']
            contig.genus = filter_df.ix[contig.contig_id]['Contig_Genus']
            contig.family = filter_df.ix[contig.contig_id]['Contig_Family']
            
        contig.calculate_signature()
        for genome in genomes:
            if contig.id == genome.id or contig.species == genome.species:
                temp_genome = deepcopy(genome)
                temp_genome.signature.subtract(contig.signature)
                temp_pseudo_par = mn.fit_nonzero_parameters([temp_genome])
                p_val = mn.log_probability(\
                    contig, temp_pseudo_par)
            else:
                p_val = mn.log_probability(\
                    contig, genome.pseudo_par)
            scores.append(\
                Score(p_val, contig, genome, contig.contig_id))

    if taxonomy_info_in_contigs:
        sys.stdout.write("p_value\tcontig_family\tcontig_genus\tcontig_species\tcontig_genome\tcompare_family\tcompare_genus\tcompare_species\tcompare_genome\tcontig_id" + os.linesep)
    for score in scores:
        sys.stdout.write(str(score) + '\n')
   
if __name__=="__main__":
    parser = ArgumentParser()
    parser.add_argument('contigs', 
                        help='specify contig file, default is stdin')
    parser.add_argument('-o', '--output', 
                        help='specify the output file.  The default is stdout')
    parser.add_argument('-t', '--taxonomy',
                        help='specify the taxonomy file.')
    parser.add_argument('-k', '--kmer_length', default=4, type=int,
                        help='specify the kmer length, default is 4')
    parser.add_argument('-d', '--directory_path',  
                        type=str, help='specify the path to where the reference genomes are located locally')
    parser.add_argument('--dir_structure', choices=['single_fasta_file','tree'], default='tree_structure',
                        help='specify the directory structure for the reference genomes')
    parser.add_argument('--no_taxonomy_info_in_contigs', action="store_true",
                        help='Add this tag if no taxonomy info is available in contigs file.')
    parser.add_argument('--filter_file',
                        help='Add filter file if not all contigs should be included')
    args = parser.parse_args()
    if args.output and args.output != '-':
        sys.stdout = open(args.output, 'w')
    if args.taxonomy:
        taxonomy_file = open(args.taxonomy, 'r')
        
    main(args.contigs, taxonomy_file, args.directory_path, args.kmer_length, args.dir_structure,not(args.no_taxonomy_info_in_contigs),filter_file=args.filter_file)

