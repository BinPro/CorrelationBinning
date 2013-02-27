#!/usr/bin/env python
from argparse import ArgumentParser
import fileinput
import numpy as np
import pandas as p
from corrbin.score import ExperimentData
from corrbin.classification import classify_bool
import matplotlib.pyplot as plt

def main(input_files, output_file):
    counter = 0
    for input_file in input_files:
        data = ExperimentData()
        data.load_data_frame(input_file)
        data.classify()
        data.calculate_roc()
    
        # Plot ROC curve
        x = data.roc_data.false_positive_rate
        y = data.roc_data.true_positive_rate
        plt.plot(x,y)
        plt.savefig(output_file + str(counter))
        counter += 1
        


if __name__=="__main__":
    parser = ArgumentParser()
    parser.add_argument('files', nargs='*', 
                        help='specify input files, default is stdin')
    parser.add_argument('-o', '--output', 
                        help='specify the base name of the output file.  The default is stdout')

    args = parser.parse_args()

    main(args.files, args.output)
