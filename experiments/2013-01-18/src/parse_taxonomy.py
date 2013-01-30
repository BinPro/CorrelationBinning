#!/usr/bin/env python
import fileinput
import sys
import os
from argparse import ArgumentParser
from collections import defaultdict

if __name__=="__main__":
    parser = ArgumentParser()
    parser.add_argument('files', nargs='*', 
        help='specify input taxa file')
    parser.add_argument('-o', '--output', 
        help='''specify the base name of the output file. 
                The default is stdout. Output will be written to
                [basename]_families.txt and [basename]_genera.txt''')
    parser.add_argument('-v', '--verbose', action='store_true',
        help='information written to stderr during execution.')
    parser.add_argument('-s', '--species', type=int, default=4,
        help='number of minimum uniq species within a genus.')
    parser.add_argument('-g', '--genera', type=int, default=4,
        help='number of minimum genera within family.')
    args = parser.parse_args()
    fh = sys.stdout
    
    families = {}
    genera = {}
    lines = []
    for line in fileinput.input(args.files):
        if not fileinput.isfirstline():
            lines.append(line)
    lines.sort(key=lambda s: s.lower())
    for line in lines:
        tmp_line = line.split('\t')
        (family,genus,species) = tmp_line[4:7]
        folder = tmp_line[0]
        if family != "NA" and genus != "NA":
            if not family in families:
                families[family] = {}
            families[family][genus] = folder
        if genus != "NA" and species !="NA":
            if not genus in genera:
                genera[genus] = {}
            genera[genus][species] = folder
    final_families = {family:family_genera for family,family_genera in families.iteritems() if len(family_genera) > args.genera}
    final_genera = {genus:genus_species for genus,genus_species in genera.iteritems() if len(genus_species) > args.species}

    if args.output and args.output != '-':
        fh = open(''.join([args.output,'_families.txt']), 'w')
    for family,genera2 in final_families.iteritems():
        print >> fh, 'group_name:\t' + family
        for genus in genera2:
            print >> fh, 'entry:\t' + genus + '\t\t' + families[family][genus]
    if args.output and args.output != '-':
        fh = open(''.join([args.output,'_genera.txt']), 'w')
    for genus,species in final_genera.iteritems():
        print >> fh, 'group_name:\t' + genus
        for specie in species:
            print >> fh, 'entry:\t' + specie + '\t\t' + genera[genus][specie]
