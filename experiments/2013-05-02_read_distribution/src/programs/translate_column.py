#!/usr/bin/env python
import fileinput
import sys
import os
from argparse import ArgumentParser
from Bio import SeqIO
import pandas as p

def main(genome_length_file,phylogeny_tsv_file,output_file,dictwise):
    genome_length_df = p.io.parsers.read_table(genome_length_file,sep='\t')
    phylogeny_df = p.io.parsers.read_table(phylogeny_tsv_file,sep='\t')
    if dictwise:
        translate_dict = {}
        fasta_names = list(phylogeny_df.fasta_name.values)
        topnames = list(phylogeny_df.topname.values)
        for i,fasta_name in enumerate(fasta_names):
            translate_dict[fasta_name] = topnames[i].strip()
        strains = genome_length_df.Strain.values
        replace_strains = []
        for strain in strains:
            try:
                replace_strains.append(translate_dict[strain])
            except KeyError:
                print >> sys.stderr, "Error, translate column is not identical to the source file column."
                print >> sys.stderr, strain + ", is not in dict"
        genome_length_df['Strain'] = p.Series(replace_strains)
        genome_length_df.to_csv(output_file,sep='\t',index=False)

if __name__=="__main__":
    parser = ArgumentParser()
    parser.add_argument('genome_length_file', 
                        help='specify genome_length file')
    parser.add_argument('-p','--phylogeny_tsv_file', 
                        help='specify phylogeny tsv file')
    parser.add_argument('-o', '--output', 
                        help='specify the output file.  The default is stdout')
    parser.add_argument('--dictwise', action='store_true',
           help='use this tag if dictwise translation should be used')

    args = parser.parse_args()
    if args.output and args.output != '-':
        sys.stdout = open(args.output, 'w')
        
    main(args.genome_length_file,args.phylogeny_tsv_file, args.output,args.dictwise)

