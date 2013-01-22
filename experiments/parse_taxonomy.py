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
        help='specify the output file.  The default is stdout')
    parser.add_argument('-v', '--verbose', action='store_true',
        help='information written to stderr during execution.')
    parser.add_argument('-s', '--species', type=int, default=4,
        help='number of minimum uniq species within a genus.')
    parser.add_argument('-g', '--genera', type=int, default=4,
        help='number of minimum genera within family.')
    args = parser.parse_args()
				      
    if args.output and args.output != '-':
        sys.stdout = open(args.output, 'w')

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
    for family,genera in final_families.iteritems():
        print "###############"
        print family
        print "###############"
        for genus in genera:
            print genus
    for genus,species in final_genera.iteritems():
        print "####################################"
        print genus
        print "####################################"
        for specie in species:
            print specie
    print len(final_families)
    print len(final_genera)
