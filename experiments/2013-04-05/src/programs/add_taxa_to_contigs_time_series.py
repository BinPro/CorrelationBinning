#!/usr/bin/env python
from argparse import ArgumentParser
import fileinput
import numpy as np
import pandas as p
from corrbin.score import ExperimentData, RocAxisFuns
from corrbin.classification import classify_bool
from corrbin.io import read_contigs_file, print_contigs_time_series, genome_info_from_parsed_taxonomy_file
from corrbin.time_series import poisson as time_series_model
import matplotlib.pyplot as plt
import sys
import probin.dna as dna


def main(contig_time_series_file, taxonomy_file,output_file):
    contig_time_series_df = p.io.parsers.read_table(contig_time_series_file,sep='\t', index_col=0)

    # Divide genomes into groups, one for each genus
    meta_genomes = genome_info_from_parsed_taxonomy_file(taxonomy_file)

    strain_meta_genome_dict = {}
    for genome_data in meta_genomes:
        strain = genome_data['species']
        strain_meta_genome_dict[strain]=genome_data

    family_dict = {}
    genus_dict = {}
    species_dict = {}
    for ix in contig_time_series_df.index:
        strain =  contig_time_series_df.ix[ix]['full_read_mappings_strain']
        family_dict[ix] = strain_meta_genome_dict[strain]['family']
        genus_dict[ix] = strain_meta_genome_dict[strain]['genus']
        species_dict[ix] = strain_meta_genome_dict[strain]['species']
        
    f_ser = p.Series(family_dict.values(),index=family_dict.keys())
    g_ser = p.Series(genus_dict.values(),index=genus_dict.keys())
    s_ser = p.Series(species_dict.values(),index=species_dict.keys())

    contig_time_series_df['Contig_Family'] = f_ser
    contig_time_series_df['Contig_Genus'] = g_ser
    contig_time_series_df['Contig_Species'] = s_ser

    contig_time_series_df.to_csv(output_file, sep='\t')
            

if __name__=="__main__":
    parser = ArgumentParser()
    parser.add_argument('file',
                        help='specify input contig time series file, default is stdin')
    parser.add_argument('-o', '--output', 
                        help='specify the base name of the output file.  The default is stdout')
    parser.add_argument('-t', '--taxonomy_file',
                        help='specify the taxonomy file')
    args = parser.parse_args()

    if args.output:
        output_file = open(args.output,'w+')
    else:
        output_file = sys.stdout
    main(args.file,open(args.taxonomy_file),output_file)

    output_file.close()
