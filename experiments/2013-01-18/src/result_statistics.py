#!/usr/bin/env python
from argparse import ArgumentParser
import sys
import pandas
import numpy as np

def main(file_name, output_base):
    df =  pandas.io.parsers.read_table(file_name, sep='\t')
    df_sts = pandas.DataFrame()
    if output_base:
        sys.stdout = open(output_base + '_species', 'w')
        
    for specie in df.contig_species.unique():
        genus = df[df.contig_species == specie].contig_genus.unique()[0]
        family = df[df.contig_species == specie].contig_family.unique()[0]

        # Self-to-self scoring:
        s_specie = df[(df.contig_species == specie) & (df.compare_species == specie)].p_value
        m_val_specie = s_specie.mean()
        std_val_specie = np.std(np.array(s_specie))
        # Within-genus scoring:
        s_genus = df[(df.contig_species == specie) & (df.compare_species != specie) & (df.compare_genus == genus)].p_value
        m_val_genus = s_genus.mean()
        std_val_genus = np.std(np.array(s_genus))
        
        # Within-family scoring:
        s_family = df[(df.contig_genus != genus) & (df.contig_species != specie) & (df.compare_genus == genus)].p_value
        m_val_family = s_family.mean()
        std_val_family = np.std(np.array(s_family))
        sys.stdout.write(specie + '\t' + str(m_val_specie) +'\t' + str(std_val_specie) + '\t' + str(m_val_genus) + '\t' + str(std_val_specie) +'\t' + str(m_val_family) +'\t' + str(std_val_family) +  '\n')

    

    # Outside-family scoring:

if __name__=="__main__":
    parser = ArgumentParser()
    parser.add_argument('file', 
                        help='specify input file, default is stdin')
    parser.add_argument('-o', '--output', 
                        help='specify the base name of the output files.  The default is stdout')
    args = parser.parse_args()

    main(args.file, args.output)
