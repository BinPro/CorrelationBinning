#!/usr/bin/env python
import fileinput
import sys
import os
from argparse import ArgumentParser
from Bio import SeqIO

def main(seq_file,read_length):
    seq_list = list(SeqIO.parse(seq_file,"fasta"))
    sys.stdout.write("\t".join(["Strain","Length"]) + "\n")
    for seq in seq_list:
        sys.stdout.write("\t".join([str(seq.id),str(len(seq))])+"\n")
    

if __name__=="__main__":
    parser = ArgumentParser()
    parser.add_argument('sequence_file', 
                        help='specify sequence file')
    parser.add_argument('-o', '--output', 
                        help='specify the output file.  The default is stdout')

    args = parser.parse_args()
    if args.output and args.output != '-':
        sys.stdout = open(args.output, 'w')
        
    main(args.sequence_file, args.output)

