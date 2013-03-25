#!/usr/bin/env python
import fileinput
import sys
import os
from copy import copy
from argparse import ArgumentParser
from probin.model.composition import multinomial as mn 
from probin.dna import DNA
from Bio import SeqIO
from corrbin.misc import all_but_index, Uniq_id, GenomeGroup
from corrbin.multinomial import Experiment
from corrbin.contig_generation import SampleSetting 
from corrbin.score import Score
from corrbin.io import read_contigs_file, genome_info_from_parsed_taxonomy_file, read_FASTA_files_no_groups

def main(contigs_file,taxonomy_file, dir_path, kmer_length):

    groups = []
    DNA.generate_kmer_hash(kmer_length)

    contigs = read_contigs_file(contigs_file)
    
    # Divide genomes into groups, one for each genus
    meta_genomes = genome_info_from_parsed_taxonomy_file(taxonomy_file)

    # Fetch sequence for each genome
    genomes = read_FASTA_files_no_groups(meta_genomes, dir_path)

    for genome in genomes:
        genome.calculate_signature()
        genome.pseudo_par = mn.fit_nonzero_parameters(genome.signature,DNA.kmer_hash_count)

    scores = []
    for contig in contigs:
        contig.calculate_signature()
        for genome in genomes:
            if contig.id == genome.id:
                temp_genome_signature = copy(genome.signature)
                temp_genome_signature.subtract(contig.signature)
                temp_pseudo_par = mn.fit_nonzero_parameters(\
                    temp_genome_signature, DNA.kmer_hash_count)
                p_val = mn.log_probability(\
                    contig.signature, temp_pseudo_par)
            else:
                p_val = mn.log_probability(\
                    contig.signature, genome.pseudo_par)
            scores.append(\
                Score(p_val, contig, genome, contig.contig_id))

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
    parser.add_argument('-d', '--directory_path', default='/home/johannes/repos/DATA/reference_genomes_ncbi', 
                        type=str, help='specify the path to where the reference genomes are located locally')
    args = parser.parse_args()
    if args.output and args.output != '-':
        sys.stdout = open(args.output, 'w')
    if args.taxonomy:
        taxonomy_file = open(args.taxonomy, 'r')
        
    main(args.contigs, taxonomy_file, args.directory_path, args.kmer_length)

