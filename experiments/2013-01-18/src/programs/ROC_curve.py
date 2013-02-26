#!/usr/bin/env python
from argparse import ArgumentParser
import fileinput
import numpy as np
import pandas as p
from corrbin.external import pyroc as roc
from corrbin.score import ExperimentData
from corrbin.classification import classify_bool
from corrbin.misc import decimal_range

def main(input_file, output_file):
    data = ExperimentData()
    data.load_data_frame(input_file)
    data.standardize()
    qs = decimal_range(-1,1,0.1)
    for q in qs:
        data.classify(q)
    


if __name__=="__main__":
    parser = ArgumentParser()
    parser.add_argument('file', 
                        help='specify input file, default is stdin')
    parser.add_argument('-o', '--output', 
                        help='specify the base name of the output file.  The default is stdout')

    args = parser.parse_args()

    main(args.file, args.output)
