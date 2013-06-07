#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu May 30 10:15:41 2013

@author: Brynjar Smari Bjarnason

Calculate statistics for known phylogenetics of contigs that have been clustered
"""
from argparse import  ArgumentParser
from corrbin import cluster_statistics as cs
from pandas import ExcelWriter
import sys

if __name__=="__main__":
    parser = ArgumentParser(description="Clustering of metagenomic contigs")
    parser.add_argument('cluster_file', 
        help='Result file from ProBin program (or each line on the form Cluster X, contig1X, contig2X,..)')
    parser.add_argument('phylo_file',
        help='Phylogenetic file for the cluster_file (each line on the form: contig_id, family, genus, species,seq_length')
    parser.add_argument('fasta_file',
        help='Fasta file to retrieve contigs lengths')
    parser.add_argument('-l', '--level', 
        default='family', type=str, choices=['family','genus','species'],
        help='Calculate statistics for which level.')
    parser.add_argument('-o', '--outfile', type=str, default=sys.stdout, 
        help='Calculate statistics for which level.')
    
    args = parser.parse_args()
    phylo_clusters_length,cm = cs.get_statistics(args.cluster_file,args.phylo_file,args.fasta_file)

    if args.outfile is not sys.stdout:
        writer = ExcelWriter(args.outfile)
        phylo_clusters_length.to_excel(writer,sheet_name="csv",header=True,index=True)
        for key,value in cm.iteritems():    
            value.to_excel(writer,sheet_name="CM {0}".format(key),header=True,index=True)
        writer.save()
    else:
        phylo_clusters_length.to_csv(args.outfile)        
    