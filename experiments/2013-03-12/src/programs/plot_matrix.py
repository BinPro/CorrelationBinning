#!/usr/bin/env python
from argparse import ArgumentParser
import fileinput
import numpy as np
import pandas as p
from corrbin.score import ExperimentData, RocAxisFuns
from corrbin.classification import classify_bool
import matplotlib.pyplot as plt

def main(input_files, output_file, title,together):
    fig = plt.figure()
    if together:
        for ix,input_file in enumerate(input_files):
            ax = plt.subplot(eval(int("13"+ ix)))
            df = p.io.parsers.read_table(input_file, sep=' & ')
            
            for contig_l in ['100','1000','10000']:
                # Plot ROC curve
                x = df.kmer_length.values
                y = df[contig_l].values
                
                label = "Contig Length: " + contig_l
                plt.plot(x,y,label=label)
    else:
        ax = plt.subplot(111)
        for input_file in input_files:
            df = p.io.parsers.read_table(input_file, sep=' & ')
            
            for contig_l in ['100','1000','10000']:
                # Plot ROC curve
                x = df.kmer_length.values
                y = df[contig_l].values
                
                label = "Contig Length: " + contig_l
                plt.plot(x,y,label=label)
    plt.ylim(0.0,1.0)
    plt.title(title)
    box = ax.get_position()
    ax.set_position([box.x0,box.y0,box.width,box.height*0.8])
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.35),
          ncol=2, fancybox=True, shadow=True)
    plt.xlabel("Kmer Length")
    plt.ylabel("Precision")
    plt.savefig(output_file)
    
        

if __name__=="__main__":
    parser = ArgumentParser()
    parser.add_argument('files', nargs='*', 
                        help='specify input files, default is stdin')
    parser.add_argument('-o', '--output', 
                        help='specify the base name of the output file.  The default is stdout')
    parser.add_argument('--title', help='Title of plot')
    parser.add_argument('--together', action="store_true", help='All matrices plotted together')
    args = parser.parse_args()

    main(args.files, args.output, args.title,args.together)
