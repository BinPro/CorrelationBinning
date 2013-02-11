#!/usr/bin/env python
from argparse import ArgumentParser
import sys
import pandas
import numpy as np

def main(file_name):
    df =  pandas.io.parsers.read_table(file_name, sep='\t')
    # Self-to-self scoring:
    df_sts = pandas.DataFrame()
    for specie in df.contig_species.unique():
        s_specie = df[(df.contig_species == specie) & (df.compare_species == specie)].p_value
        m_val = s_specie.mean()
        std_val = np.std(np.array(s_specie))
        sys.stdout.write(specie + '\t' + str(m_val) +'\t' + str(std_val) + '\n')

    # Within-genus scoring:
    
    # Within-family scoring:

    # Outside-family scoring:

if __name__=="__main__":
    parser = ArgumentParser()
    parser.add_argument('file', 
                        help='specify input file, default is stdin')
    parser.add_argument('-o', '--output', 
                        help='specify the output file.  The default is stdout')
    args = parser.parse_args()
    if args.output and args.output != '-':
        sys.stdout = open(args.output, 'w')
        
    main(args.file)
