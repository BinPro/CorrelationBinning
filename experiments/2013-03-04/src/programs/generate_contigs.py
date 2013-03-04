#!/usr/bin/env python
import fileinput
import sys
import os
from argparse import ArgumentParser
from probin.dna import DNA
from Bio import SeqIO
from corrbin.misc import all_but_index, Uniq_id
from corrbin.contig_generation import SampleSetting, sample_contig

def main(open_name_file, dir_path, x_set):

    DNA.generate_kmer_hash(kmer_length)

    groups = read_parsed_taxonomy_file(open_name_file)

    # Read the FASTA files for each genome
    
    # Each genome in a group is a bin, fit parameters to all bins
    os.chdir(dir_path)
    for group in groups:
        for genome_data in group.genome_data:
            dir_name = genome_data['file_name']
            fasta_files = os.listdir(dir_name)
            for fasta_file in fasta_files:
                genome_file = open(dir_name + '/' + fasta_file)
                identifier = genome_file.readline()
                # Only use non-plasmid genomes
                # Some bacterial genomes contain more than 1 chromosonme,  
                # but assumed not more than 2
                if identifier.find('plasmid') == -1 and identifier.find('chromosome 2') == -1:
                    genome_file.close() #Close and reopen the same file
                    genome_file = open(dir_name + '/' + fasta_file)
                    genome_seq = list(SeqIO.parse(genome_file, "fasta"))
                    if len(genome_seq) > 1:
                        sys.stderr.write("Warning! The file " + fasta_file + " in directory " + dir_name + " contained more than one sequence, ignoring all but the first!" + os.linesep)
                    genome = DNA(id = dir_name, seq= str(genome_seq[0].seq))
                    genome.genus = genome_data['genus']
                    genome.species = genome_data['species']
                    genome.family = genome_data['family']
                    group.genomes.append(genome)
                genome_file.close()

    # For each bin, generate a number of contigs, 
    # re-calculate parameters for that bin without contig-section.
    # Further score this contig against all bins, keep within-group
    # scores separate from outside-group scores.
    all_scores = []
    id_generator = Uniq_id(1000)
    for group_index in range(len(groups)):
        group = groups[group_index]
        rest_groups = all_but_index(groups, group_index)
        test = Test(x_set, group, rest_groups, id_generator)
        group_scores = test.execute()
        
        all_scores.append(group_scores)
    sys.stdout.write("p_value\tcontig_family\tcontig_genus\tcontig_species\tcontig_genome\tcompare_family\tcompare_genus\tcompare_species\tcompare_genome\tcontig_id" + os.linesep)
    for group_scores in all_scores:
        for genome_scores in group_scores:
            for score in genome_scores:
                sys.stdout.write(str(score) + '\n')
 

if __name__=="__main__":
    parser = ArgumentParser()
    parser.add_argument('file',
                        help='specify input file, containing parsed genome data, default is stdin')
    parser.add_argument('-o', '--output', 
                        help='specify the output file.  The default is stdout')
    parser.add_argument('-d', '--directory_path', default='/home/johannes/repos/DATA/reference_genomes_ncbi', 
                        type=str, help='specify the path to where the reference genomes are located')
    parser.add_argument('-c', '--no_contigs', default=100, type=int,
                        help='Specify the number of contigs to be sampled from each group. This may be only approximate due to what priority is chosen') 
    parser.add_argument('-p', '--priority', default="genomes",
                        type=str, help='specify the prioritized way of sampling contigs. Specify "groups" to make sure each group is sampled exactly the number of times specified by no_contigs, distributed randomly over the genomes present in each group, or specify "genomes" to make sure each genome within a certain group contributes with exactly the same number of contigs. This setting makes -c setting only approximate.')
    parser.add_argument('--contig_min_length', default=1000, type=int, help='Specify the minimum length for contigs')
    parser.add_argument('--contig_max_length', default=1000, type=int, help='Specify the maximum length for contigs')
    parser.add_argument('--debug_mode', action='store_true', help='In debug mode, all contigs will start at the first nucleotide, making testing possible.')
    args = parser.parse_args()
    if args.output and args.output != '-':
        sys.stdout = open(args.output, 'w')
        
    name_file_handle = fileinput.input(args.file)
    sample_setting = SampleSetting(args.priority, args.no_contigs, args.contig_min_length, args.contig_max_length, args.debug_mode)
    main(name_file_handle, args.directory_path, args.kmer_length, sample_setting)
    name_file_handle.close()
