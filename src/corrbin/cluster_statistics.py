#!/usr/bin/env python
"""
Created on Thu Mar 28 10:07:24 2013
Calculate recall and precision for clustering algorithms. Clustering 
@author: Brynjar Smari Bjarnason
"""
import sys
import os
import numpy as np
from pandas import DataFrame, Series, pivot_table, read_csv


from argparse import ArgumentParser

def get_phylo_from_file(phylo_file):
    """
    phylo_file: File in path is a comma seperated file with columns: id, family, genus, specie
                family      genus       species    
    contig1     family1     genus1      species1
    contig2     family1     genus2      species2
    contig1X    family2     genus3      species3
    contig2X    family2     genus3      species4
    ...
    """
    df = read_csv(phylo_file,index_col=0)
    return df

def get_clusters_from_file(cluster_file):
    """
    cluster_file: File with the output from ProBin execution. Basically, 
    the clusters should be represented as follows:
    Cluster 0, contig1, contig2, ...
    Cluster 1, contig1X, contig2X, ...
    ...
    result: Padas data frame with column cluster and the contig ids as indexes.
                cluster    
    contig1     cluster0
    contig2     cluster0
    contig1X    cluster1
    contig2X    cluster1
    ...
    """
    index=[]
    values=[]
    with open(cluster_file) as fh:
        for clustering in fh.readlines():
            if not clustering.startswith("#"):
                s = clustering.strip().split(",")
                cluster = s[0]
                contigs = s[1:]
                for contig in contigs:
                    index.append(contig)
                    values.append(cluster)
    clustering = DataFrame(Series(values,index=index),columns=["cluster"])
                
    return clustering

def get_statistics(cluster_file, phylo_file,level):
    clusters = get_clusters_from_file(cluster_file)
    phylo = get_phylo_from_file(phylo_file)
    if len(clusters) != len(phylo):
        print >> sys.stderr, "not equally many contigs in clustering({0}) and phylo({1}).".format(len(clusters),len(phylo))
        sys.exit(-1)
    #we can do this join on the dataframes since the indexes are the same contigs!
    phylo_clusters = phylo.join(clusters)
    
    return phylo_clusters


#    _get_phylo(contigs)    
#    cm = confusion_matrix(contigs,clusters)
#    rc = recall(contigs,clusters)
#    pr = precision(contigs,clusters)
#    for c,r,p in izip(cm,rc,pr):
#        print>>output,"#{0}".format("="*20)
#        c.to_csv(output)
#        r.recall.to_csv(output)
#        p.precision.to_csv(output)

def confusion_matrix(contigs,clustering,level="family"):
    levels = ["family","genus","species"]
    df = _create_dataframe(contigs,clustering)
    cm = [pivot_table(df,rows=levels[:i+1],cols=["cluster"],aggfunc=np.sum) for i in xrange(len(levels))]
    return cm

def recall(contigs,clustering):
    confusion_matrixes = confusion_matrix(contigs,clustering)
    recalls = [matrix.div(matrix.sum(axis=1),axis=0) for matrix in confusion_matrixes]
    for matrix in recalls:
        matrix["recall"] = matrix.max(axis=1)
    return recalls

def precision(contigs,clustering):
    confusion_matrixes = confusion_matrix(contigs,clustering)
    precisions = [matrix.div(matrix.sum()).T for matrix in confusion_matrixes]
    for precision in precisions:
        precision["precision"] = precision.max(axis=1)
    return precisions