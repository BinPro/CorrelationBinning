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
import re
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from collections import defaultdict

if __name__=="__main__":
    parser = ArgumentParser(description="Clustering of metagenomic contigs")
    parser.add_argument('cluster_files', nargs="+", 
        help='Result file from ProBin program (or each line on the form Cluster X, contig1X, contig2X,..)')
    parser.add_argument('phylo_file',
        help='Phylogenetic file for the cluster_file (each line on the form: contig_id, family, genus, species,seq_length')
    parser.add_argument('-l', '--level', 
        default='family', type=str, choices=['family','genus','species'],
        help='Calculate statistics for which level.')
    parser.add_argument('-o', '--outfile', action='store_true', default=False,
        help='Calculate statistics for which level.')
    args = parser.parse_args()
    results = [(cs.get_statistics(cluster_file,args.phylo_file)+(cluster_file,)) for cluster_file in args.cluster_files]
    for (phylo_clusters_length,cm,cluster_file) in results:
        if args.outfile:
            writer = ExcelWriter("{0}.{1}".format(args.cluster_files[0],"xlsx"))
            phylo_clusters_length.to_excel(writer,sheet_name="csv",header=True,index=True)
            for key,value in cm.iteritems():    
                value[0].to_excel(writer,sheet_name="CM {0}".format(key),header=True,index=True)
        writer.save()
    else:
        phylo_clusters_length.to_csv("{0}.{1}".format(args.cluster_files[0],"csv"))
#==============================================================================
#    cluster_level = defaultdict(dict)
#    for phylo_clusters_length, cm, cluster_file in results:
#        re_clusters = re.compile("c(\d+)")
#        re_kmer = re.compile("k(\d+)")
#        re_alg = re.compile("(em|kmeans)")
#        re_contig = re.compile("100_(100\d*)")
#
#        nr_clusters = re.search(re_clusters,cluster_file).groups()[0]
#        kmers = re.search(re_kmer,cluster_file).groups()[0]
#        alg = re.search(re_alg,cluster_file).groups()[0]
#        contig = re.search(re_contig,cluster_file).groups()[0]
#        
#        clust_alg = []
#        clust_alg.append((13,cm["family"][1],cm["family"][2]))
#        clust_alg.append((55,cm["genus"][1],cm["genus"][2]))
#        clust_alg.append((184,cm["species"][1],cm["species"][2]))
#        cluster_level[nr_clusters][alg] = clust_alg
#        
#    for k,algs in cluster_level.iteritems():
#        em_values = algs["em"]
#        kmeans_values = algs["kmeans"]
#        fig = plt.figure()
#        
#        ax = fig.add_subplot(111,title="n=%s" % k)
#        plt.ylabel("ratio")
#        plt.xlabel("nr_clusters")
#        ax.plot([x[0] for x in em_values], [x[1] for x in em_values], "bo-",label='em_prec')
#        ax.plot([x[0] for x in em_values], [x[2] for x in em_values],"rx-",label= 'em_rec')
#        ax.plot([x[0] for x in kmeans_values], [x[1] for x in kmeans_values], "go-",label="kmeans_prec")
#        ax.plot([x[0] for x in kmeans_values], [x[2] for x in kmeans_values],"yx-",label="kmeans_rec")
#        ax.legend(bbox_to_anchor=(0., 1.10, 1., .102), loc=3,ncol=2, mode="expand", borderaxespad=0.)
#        plt.ylim(0,1)
#        plt.savefig("/tmp/contig_{2}_k_{0}_c_{1}.png".format(kmers,k,contig),bbox_inches='tight')
#==============================================================================
