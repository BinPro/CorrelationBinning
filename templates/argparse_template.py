#!/usr/bin/env python
import fileinput
import sys
import os
from argparse import ArgumentParser

from Bio import SeqIO

def main():
    sys.stdout.write("Hello world! \n")

if __name__=="__main__":
    parser = ArgumentParser()
    parser.add_argument('files', nargs='*', 
        help='specify input files, default is stdin')
    parser.add_argument('-o', '--output', 
        help='specify the output file.  The default is stdout')
    parser.add_argument('-v', '--verbose', action='store_true',
        help='information written to stderr during execution.')
    parser.add_argument('-n', '--number', nargs='?', default=100, type=int,
        help='specify the number of contigs SOMETHING, default is 100')
    args = parser.parse_args()
    if args.output and args.output != '-':
        sys.stdout = open(args.output, 'w')

    handle = fileinput.input(args.files)
    genomes = list(SeqIO.parse(handle,"fasta"))
    handle.close()
    if args.verbose:
        sys.stderr.write("Number of genomes read: %i %s" % (len(genomes),os.linesep))
    
    main()
