#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu May 30 10:15:41 2013

@author: Brynjar Smari Bjarnason

Calculate statistics for known phylogenetics of contigs that have been clustered
"""
from argparse import FileType, ArgumentParser
from corrbin import cluster_statistics as cs
import sys

if __name__=="__main__":
    parser = ArgumentParser(description="Clustering of metagenomic contigs")
    parser.add_argument('cluster_file', 
        help='Result file from ProBin program (or each line on the form Cluster X, contig1X, contig2X,..)')
    parser.add_argument('phylo_file',
        help='Phylogenetic file for the cluster_file (each line on the form: contig_id, family, genus, species')
    parser.add_argument('-l', '--level', 
        default='family', type=str, choices=['family','genus','species'],
        help='Calculate statistics for which level.')
    parser.add_argument('-o', '--outfile', type=FileType('w'), default=sys.stdout, 
        help='Calculate statistics for which level.')
    
    args = parser.parse_args()
    phylo_clusters = cs.get_statistics(args.cluster_file,args.phylo_file,args.level)
    phylo_clusters.to_csv(args.outfile)