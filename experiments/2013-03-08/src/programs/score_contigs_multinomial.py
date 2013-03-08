#!/usr/bin/env python
import fileinput
import sys
import os
from argparse import ArgumentParser
from probin.model.composition import multinomial as mn 
from probin.dna import DNA
from Bio import SeqIO
from corrbin.misc import all_but_index, Uniq_id, GenomeGroup
from corrbin.multinomial import Experiment
from corrbin.contig_generation import SampleSetting
from corrbin.score import read_contigs_file

def main(contigs_file,taxonomy_file, dir_path, kmer_length):

    groups = []
    DNA.generate_kmer_hash(kmer_length)

    contigs = read_contigs_file(contigs_file)
    
    # Divide genomes into groups, one for each genus
    groups = read_parsed_taxonomy_file(taxonomy_file)

    # Fetch sequence for each genome
    genomes = read_FASTA_files_no_groups(genome_names, dir_path)

    for genome in genomes:
        genome.pseudo_par = genome.fit_nonzero_parameters(genome.sig,DNA.kmer_hash_count)

    for contig in contigs:
        for genome in genomes:
            p_val = mn.log_probability(contig.signature, genome.pseudo_par)
            scores.append(Score(p_val, contig, genome, contig.contig_id))

    sys.stdout.write("p_value\tcontig_family\tcontig_genus\tcontig_species\tcontig_genome\tcompare_family\tcompare_genus\tcompare_species\tcompare_genome\tcontig_id" + os.linesep)
    for score in scores:
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
    ex_setting = SampleSetting(args.priority, args.no_contigs, args.contig_min_length, args.contig_max_length, args.debug_mode)
    main(name_file_handle, args.directory_path, args.kmer_length, ex_setting)
    name_file_handle.close()
