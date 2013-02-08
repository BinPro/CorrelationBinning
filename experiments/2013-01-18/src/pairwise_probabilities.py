#!/usr/bin/env python
import fileinput
import sys
import os
from argparse import ArgumentParser
from probin.model.composition import multinomial as mn 
from probin.dna import DNA
from Bio import SeqIO
from helpers import all_but_index, GenomeGroup, ExperimentSetting, Test

def main(open_name_file, dir_path, kmer_length, x_set):

    groups = []
    DNA.generate_kmer_hash(kmer_length)
    # Read the file with all names, divide them into groups
    for line in open_name_file:
        if line[0:12] == 'family_name:':
            family = line.split('\t')[1].strip()
        elif line[0:11] == 'genus_name:':
            genus = line.split('\t')[1].strip()
            new_group = GenomeGroup(genus)
            new_group.family = family
            groups.append(new_group)
        elif line[0:6] == 'entry:':
            genome_name = line.split('\t')[2].strip()
            genome_species = line.split('\t')[1].strip()
            meta_genome = {'id': genome_name,
                           'species': genome_species,
                           'genus': genus,
                           'family': family,
                           'file_name': genome_name
                          }
            groups[-1].genome_data.append(meta_genome)

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
    for group_index in range(len(groups)):
        group = groups[group_index]
        rest_groups = all_but_index(groups, group_index)
        test = Test(x_set, group, rest_groups)
        group_scores = test.execute()
        
        all_scores.append(group_scores)
    for group_scores in all_scores:
        for genome_scores in group_scores:
            for score in genome_scores:
                sys.stdout.write(str(score) + '\n')
 

if __name__=="__main__":
    parser = ArgumentParser()
    parser.add_argument('files', nargs='*', 
                        help='specify input files, default is stdin')
    parser.add_argument('-o', '--output', 
                        help='specify the output file.  The default is stdout')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='information written to stderr during execution.')
    parser.add_argument('-m', '--model', default="multinomial", type=str,
                        help='specify the model to use for calculating the probabilities, default is multinomial')
    parser.add_argument('-k', '--kmer_length', default=4, type=int,
                        help='specify the kmer length, default is 4')
    parser.add_argument('-d', '--directory_path', default='/home/johannes/repos/DATA/reference_genomes_ncbi', 
                        type=str, help='specify the path to where the reference genomes are located locally')
    parser.add_argument('-c', '--no_contigs', default=100, type=int,
                        help='Specify the number of contigs to be sampled from each group. This may be only approximate due to what priority is chosen') 
    parser.add_argument('-p', '--priority', default="genomes",
                        type=str, help='specify the prioritized way of sampling contigs. Specify "groups" to make sure each group is sampled exactly the number of times specified by no_contigs, distributed randomly over the genomes present in each group, or specify "genomes" to make sure each genome within a certain group contributes with exactly the same number of contigs.')
    parser.add_argument('--contig_min_length', default=1000, type=int, help='Specify the minimum length for contigs')
    parser.add_argument('--contig_max_length', default=1000, type=int, help='Specify the maximum length for contigs')
    parser.add_argument('--debug_mode', action='store_true', help='In debug mode, all contigs will start at the first nucleotide, making testing possible.')
    args = parser.parse_args()
    if args.output and args.output != '-':
        sys.stdout = open(args.output, 'w')
        
    name_file_handle = fileinput.input(args.files)
    if args.verbose:
        sys.stderr.write("Number of genomes read: %i %s" % (len(genomes),os.linesep))
    mode = "refit"
    ex_setting = ExperimentSetting(args.priority, mode, args.no_contigs, args.contig_min_length, args.contig_max_length, args.debug_mode)
    main(name_file_handle, args.directory_path, args.kmer_length, ex_setting)
    name_file_handle.close()
