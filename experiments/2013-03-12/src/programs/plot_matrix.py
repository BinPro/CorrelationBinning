#!/usr/bin/env python
from argparse import ArgumentParser
import fileinput
import numpy as np
import pandas as p
from corrbin.score import ExperimentData, RocAxisFuns
from corrbin.classification import classify_bool
import matplotlib.pyplot as plt
import sys

def main(input_files, output_file, title,together,mock_multiple_samples):
    if together:
        fig = plt.figure(figsize=(15,6))
        lines =[]
        labels = []
        for ix,input_file in enumerate(input_files):
            if ix == 0:
                ax0 = plt.subplot(1,3,ix+1)
                ax = ax0
            else:
                ax = plt.subplot(1,3,ix+1,sharey=ax0)
            df = p.io.parsers.read_table(input_file, sep=' & ')

            linestyles = [':','--','-.']
            for contig_nr,contig_l in enumerate(['100','1000','10000']):
                x = df.kmer_length.values
                y = df[contig_l].values
                
                label = "Contig Length: " + contig_l
                labels.append(label)
                line = plt.plot(x,y)
                lines.append(line)
                plt.setp(line,linestyle=linestyles[contig_nr],linewidth=2.0, label=label)
            if ix == 2:
                plt.setp(ax.get_yticklabels(),visible=False)
                plt.setp(ax.get_xticklabels(),fontsize=14)
            elif ix == 1:
                plt.setp(ax.get_yticklabels(),visible=False)
                plt.setp(ax.get_xticklabels(),fontsize=14)
                plt.xlabel("Kmer Length",fontsize=16)           
            else:
                plt.setp(ax.get_yticklabels(),fontsize=14)
                plt.setp(ax.get_xticklabels(),fontsize=14)
                plt.ylabel("Precision",fontsize=16)


            plt.ylim(0.0,1.0)
            box = ax.get_position()
            ax.set_position([box.x0,box.y0,box.width,box.height*0.8])

    elif mock_multiple_samples:
        fig = plt.figure()
        ax = plt.subplot(111)
        df = p.io.parsers.read_table(input_file,sep=' & ')
        x = df.nr_samples.values
        y = df.values
        plt.plot(x,y)
    else:
        fig = plt.figure()
        ax = plt.subplot(111)
        for input_file in input_files:
            df = p.io.parsers.read_table(input_file, sep=' & ')
            
            for contig_l in ['100','1000','10000']:
                # Plot ROC curve
                x = df.kmer_length.values
                y = df[contig_l].values
                
                label = "Contig Length: " + contig_l
                plt.plot(x,y,label=label)

    plt.legend(loc='upper center', bbox_to_anchor=(-0.7, 1.25),
               ncol=2, fancybox=True, shadow=True)
    plt.savefig(output_file)
    
        

if __name__=="__main__":
    parser = ArgumentParser()
    parser.add_argument('files', nargs='*', 
                        help='specify input files, default is stdin')
    parser.add_argument('-o', '--output', 
                        help='specify the base name of the output file.  The default is stdout')
    parser.add_argument('--title', help='Title of plot')
    parser.add_argument('--together', action="store_true", help='All matrices plotted together')
    parser.add_argument('--mock_multiple', action="store_true", help='Plot multimple samples from mock data')
    args = parser.parse_args()

    main(args.files, args.output, args.title,args.together, args.mock_multiple)
