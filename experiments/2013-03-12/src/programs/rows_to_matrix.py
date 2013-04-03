#!/usr/bin/env python
from argparse import ArgumentParser
import fileinput
import numpy as np
import pandas as p
from corrbin.score import ExperimentData, RocAxisFuns
from corrbin.classification import classify_bool
import matplotlib.pyplot as plt
from glob import glob
import sys

def main(input_files):
    fig = plt.figure()
    ax = plt.subplot(111)

    file_names = {}
    for input_file in input_files:
        level = input_file.split("_")[-1]
        kmer_len = int(input_file.split("_")[-2])
        open_input = open(input_file)
        content = open_input.read()
        content += '\n'
        if level in file_names:
            if kmer_len in file_names[level]:
                pass
            else:
                file_names[level][kmer_len] = content
        else:
            file_names[level] = {kmer_len: content}

    for level in file_names.keys():
        output_file = open(level + "_matrix",'w+')
        output_file.write(" & ".join(["kmer_length","100","1000","10000"])+"\n")
        for kmer_len, content in iter(sorted(file_names[level].iteritems())):
            output_file.write(str(kmer_len) +" & " + content)


    

if __name__=="__main__":
    parser = ArgumentParser()
    parser.add_argument('files', nargs='*', 
                        help='specify input files, default is stdin')

    args = parser.parse_args()

    main(args.files)
