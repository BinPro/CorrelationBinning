#!/usr/bin/env python
import fileinput
import sys
import os
from argparse import ArgumentParser
from collections import defaultdict
from copy import deepcopy

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
    parser.add_argument('--file_name_column', default='genome',
        help='specify the name of the column containing the file_name, default genome')
    parser.add_argument('--backup_file_name_column',
        default='genome', help='specify backup file name column if standard == "N/A"')      
    parser.add_argument('--species_column', default='species',
        help='specify the name of the column containing the species name, default species')
    
    args = parser.parse_args()
    fh = sys.stdout
    
    families = {}
    genera = {}
    lines = []
    for line in fileinput.input(args.files):
        if not fileinput.isfirstline():
            lines.append(line)
        else:
            keys = line.strip().split('\t')
    lines.sort(key=lambda s: s.lower())
    for line in lines:
        tmp_line = line.strip().split('\t')
        tmp_hash = {}
        i=0
        for key in keys:
            tmp_hash[key] = tmp_line[i]
            i +=1

        family = tmp_hash['family']
        genus = tmp_hash['genus']
        species = tmp_hash['species']
        if family != "NA" and genus != "NA" and species !="NA":
            if not family in families:
                families[family] = {}
                families[family][genus] = {}
            elif not genus in families[family]:
                families[family][genus] = {}
            families[family][genus][species] = tmp_hash
    # remove genera with too few species
    mid_families = deepcopy(families)
    for family, genera_hash in families.iteritems():
        for genera, species_hash in genera_hash.iteritems():
            if len(species_hash) < args.species:
                mid_families[family].pop(genera)
    final_families = {family:family_genera for family,family_genera in mid_families.iteritems() if len(family_genera) > args.genera}

    if args.output and args.output != '-':
        fh = open(''.join([args.output,'_complete.txt']), 'w')
    for family,genera2 in final_families.iteritems():
        print >> fh, 'family_name:\t' + family
        for genus in genera2:
            print >> fh, 'genus_name:\t' + genus
            for specie,tmp_hash in families[family][genus].iteritems():
                file_name = families[family][genus][specie][args.file_name_column]
                if file_name == "N/A":
                    file_name = families[family][genus][specie][args.backup_file_name_column]

                print >> fh, 'entry:\t' + specie + '\t' + file_name                    
