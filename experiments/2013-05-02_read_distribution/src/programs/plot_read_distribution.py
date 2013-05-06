#!/usr/bin/env python
import fileinput
import sys
import os
import numpy as np
from argparse import ArgumentParser
from probin.model.coverage import binomial as model
from probin.dna import DNA
import matplotlib.pyplot as plt
import pandas as p

def main(read_mapping_file_name,read_length,output_file):
    df_rm = p.io.parsers.read_table(read_mapping_file_name,sep='\t',index_col=0)
    cl = df_rm[df_rm.full_read_mappings != 'N/A'].contig_length
    frm= df_rm[df_rm.full_read_mappings != 'N/A'].full_read_mappings.astype(float)
    df_rm['coverage'] = frm*read_length/cl
    strains = df_rm.full_read_mappings_strain.unique()
    first_strain_df = df_rm
    x = np.log(first_strain_df.contig_length.values)
    y = first_strain_df[np.isfinite(first_strain_df.coverage) & (first_strain_df.coverage < 15)].coverage.values
    
    sys.stderr.write(str(np.mean(y))+'\n')
#    z = np.random.exponential(scale=np.mean(y),size=len(y))
    z = np.random.poisson(lam=np.mean(y)*0.1,size=len(y))/0.1
    sys.stderr.write(str(len(y))+"SPACE"+str(np.max(y))+'\n')
    sys.stderr.write(str(len(z))+"SPACE"+str(np.max(z))+'\n')
    plt.hist(y,bins=200,color='b',label="Coverage")
    plt.hist(z,bins=200,color='r',label="Poisson",alpha=0.5,range=(0,15))
    plt.savefig(output_file)    
    

if __name__=="__main__":
    parser = ArgumentParser()
    parser.add_argument('read_mapping_file', 
                        help='specify read mapping file with start positions, default is stdin')
    parser.add_argument('-o', '--output', 
                        help='specify the output file.  The default is stdout')
    parser.add_argument('-l', '--read_length', default=100, type=int,
                        help='specify the read length, default is 100')
    args = parser.parse_args()
        
    main(args.read_mapping_file, args.read_length, args.output)

