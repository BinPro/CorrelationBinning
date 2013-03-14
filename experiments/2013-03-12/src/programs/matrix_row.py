#!/usr/bin/env python
from argparse import ArgumentParser
import fileinput
import sys
import numpy as np
import pandas as p
from corrbin.score import ExperimentData, RocAxisFuns, RocAxisFun
from corrbin.classification import classify_bool
import matplotlib.pyplot as plt

def main(input_files, level, fun_name):
    fun = RocAxisFun(fun_name).fun
    m_row = []

    for input_file in input_files:
        data = ExperimentData(fun)
        data.load_data_frame(input_file)
        data.standardize()
        data.classify(level)
        df = data.classification[level]
        real_classif = [row[1] for row in df.values]
        m_row.append(fun(real_classif,[]))

    formated_matrix_row = " & ".join([str(n) for n in m_row])
    sys.stdout.write(formated_matrix_row)


if __name__=="__main__":
    parser = ArgumentParser()
    parser.add_argument('files', nargs='*', 
                        help='specify input files, default is stdin')
    parser.add_argument('-o', '--output', 
                        help='specify the base name of the output file.  The default is stdout')
    parser.add_argument('--level',
                        help='specify levels to investigate, e.g. family')
    parser.add_argument('-f','--function', choices=RocAxisFuns.FUNS,
                        help='The function statistic')

    args = parser.parse_args()

    default_stdout = sys.stdout
    if args.output and args.output != '-':
        sys.stdout = open(args.output,'w')

    main(args.files, args.level, args.function)

    sys.stdout = default_stdout
