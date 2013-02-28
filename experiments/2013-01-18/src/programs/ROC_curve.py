#!/usr/bin/env python
from argparse import ArgumentParser
import fileinput
import numpy as np
import pandas as p
from corrbin.score import ExperimentData
from corrbin.classification import classify_bool
import matplotlib.pyplot as plt

def main(input_files, output_file, levels):
    fig = plt.figure(1)
    for input_file in input_files:
        data = ExperimentData()
        data.load_data_frame(input_file)
        for level in levels:
            data.classify(level)
        data.calculate_roc()
    
        # Plot ROC curve
        for level in levels:
            x = data.roc_data[level].false_positive_rate
            y = data.roc_data[level].true_positive_rate
        
            file_name = input_file.split('/')[-1]
            plt.plot(x,y,label=file_name + ", " + level)
    
    plt.legend()
    plt.savefig(output_file)
    
        


if __name__=="__main__":
    parser = ArgumentParser()
    parser.add_argument('files', nargs='*', 
                        help='specify input files, default is stdin')
    parser.add_argument('-o', '--output', 
                        help='specify the base name of the output file.  The default is stdout')
    parser.add_argument('--levels', nargs='*',
                        help='specify levels to investigate, e.g. family')

    args = parser.parse_args()

    main(args.files, args.output, args.levels)
