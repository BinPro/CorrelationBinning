#!/usr/bin/env python
from argparse import ArgumentParser
import fileinput
import sys
import numpy as np
import matplotlib.pyplot as plt
import pandas as p


def main(input_file, output_file, min_score, max_score):
    hist_range = (min_score,max_score)
    df =  p.io.parsers.read_table(input_file, sep='\t')
    df_self_to_self = df[df.contig_species == df.compare_species]
    df_within_genus = df[(df.contig_species != df.compare_species) & (df.contig_genus == df.compare_genus)]
    df_within_family = df[(df.contig_genus != df.compare_genus) & (df.contig_family == df.compare_family)]
    fig = plt.figure(1)
    plt.subplot(311,title="Self to self")
    df_self_to_self = df_self_to_self[np.isfinite(df_self_to_self.p_value)]
    plt.hist(df_self_to_self.p_value, 30, hist_range)
    plt.subplot(312, title="Within Genus")
    df_within_genus = df_within_genus[np.isfinite(df_within_genus.p_value)]
    plt.hist(df_within_genus.p_value,30, hist_range)
    plt.subplot(313, title="Within Family")
    df_within_family = df_within_family[np.isfinite(df_within_family.p_value)]
    plt.hist(df_within_family.p_value,30, hist_range)
    fig.tight_layout()

    plt.savefig(output_file)


if __name__=="__main__":
    parser = ArgumentParser()
    parser.add_argument('file', 
                        help='specify input file, default is stdin')
    parser.add_argument('-o', '--output', 
                        help='specify the base name of the output file.  The default is stdout')
    parser.add_argument('--min_score', type=int,
                        help='specify the left limit of the histograms')
    parser.add_argument('--max_score', type=int,
                        help='specify the right limit of the histograms')
    args = parser.parse_args()

    main(args.file, args.output, args.min_score, args.max_score)
