#!/usr/bin/env python
import fileinput
import sys
import os
from argparse import ArgumentParser
from probin.dna import DNA
from Bio import SeqIO
from corrbin.misc import all_but_index, Uniq_id, GenomeGroup
from corrbin.io import read_parsed_taxonomy_file, \
    read_FASTA_files, print_parts

def main(open_name_file, dir_path, n, l):

    DNA.generate_kmer_hash(1)

    groups = read_parsed_taxonomy_file(open_name_file)

    # Read in the FASTA files for each genome
    read_FASTA_files(groups,dir_path)

    # For each bin, generate a number of contigs, 
    all_scores = []
    id_generator = Uniq_id(1000)
    for group_index,group in enumerate(groups):
        for genome in group.genomes:
            parts = genome.split_seq_random(l,n)
            print_parts(parts,sys.stdout, id_generator, genome)


if __name__=="__main__":
    parser = ArgumentParser()
    parser.add_argument('file',
                        help='specify input file, containing parsed genome data, default is stdin')
    parser.add_argument('-o', '--output', 
                        help='specify the output file.  The default is stdout')
    parser.add_argument('-d', '--directory_path', default='/home/johannes/repos/DATA/reference_genomes_ncbi', 
                        type=str, help='specify the path to where the reference genomes are located')
    parser.add_argument('-n', '--number_of_parts', default=100, type=int,
                        help='Specify the number of parts to be sampled from each genome.') 
    parser.add_argument('--part_length', default=1000, type=int, help='Specify the length for the parts')
    args = parser.parse_args()
    if args.output and args.output != '-':
        sys.stdout = open(args.output, 'w')
        
    taxonomy_file_handle = fileinput.input(args.file)
    main(taxonomy_file_handle, args.directory_path, args.number_of_parts, args.part_length)
    taxonomy_file_handle.close()
    sys.stdout.close()
