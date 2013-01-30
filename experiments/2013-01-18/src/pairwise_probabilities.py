#!/usr/bin/env python
import fileinput
import sys
import os
from argparse import ArgumentParser

from Bio import SeqIO

def main(open_name_file, dir_path):
    groups = []
    # Read the file with all names, divide them into groups
    for line in open_name_file:
        if line[0:11] == 'group_name:':
            groups.append((line.split('\t')[1].strip(), []))
        elif line[0:6] == 'entry:':
            groups[-1][-1].append(line.split('\t')[3].strip())

    # Each genome in a group is a bin, fit parameters to all bins
    os.chdir(dir_path)
    for group in groups:
        for dir_name in group[-1]:
            print dir_name
            fasta_files = os.listdir(dir_name)
            for fasta_file in fasta_files:
                genome_file = open(dir_name + '/' + fasta_file)
                identifier = genome_file.readline()
                # Only use non-plasmid genomes
                if identifier.find('plasmid') == -1 and identifier.find('chromosome 2') == -1:
                    #Does not work: genome = list(SeqIO.parse(genome_file, "fasta"))
                    print fasta_file
                genome_file.close()

            
            

    # For each bin, generate a contig, re-calculate parameters for
    # that bin without contig-section.

    # Score this contig against all bins, keep within-group
    # scores separate from outside-group scores.

 

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
    args = parser.parse_args()
    if args.output and args.output != '-':
        sys.stdout = open(args.output, 'w')

    name_file_handle = fileinput.input(args.files)
    if args.verbose:
        sys.stderr.write("Number of genomes read: %i %s" % (len(genomes),os.linesep))
        
    main(name_file_handle, args.directory_path)
    name_file_handle.close()
