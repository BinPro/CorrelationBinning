"""
Created on Thu Mar 28 10:07:24 2013
Calculate recall and precision for clustering algorithms. Clustering 
@author: Brynjar Smari Bjarnason
"""
import pandas as pd
from pandas import DataFrame, Series, pivot_table
import numpy as np
from itertools import izip
from argparse import ArgumentParser

def get_phylo_from_file(phylo_file):
    """
    File in path is a comma seperated file with columns: contig_id, family, genus, species
    """
    return pd.read_csv(phylo_file)

def get_clusters_from_file(cluster_file):
    """
    File with the output from ProBin execution. Basically, the clusters should be represented as follows:
    Cluster 0, contig1, contig2, ...
    Cluster 1, contig1X, contig2X, ...
    ...
    """
    return pd.read_csv(cluster_file)

def get_statistics(cluster_file, phylo_file,level):
    clusters = get_clusters_from_file(cluster_file)
    phylo = get_phylo_from_file(phylo_file)

    print phylo
    print clusters


#    _get_phylo(contigs)    
#    cm = confusion_matrix(contigs,clusters)
#    rc = recall(contigs,clusters)
#    pr = precision(contigs,clusters)
#    for c,r,p in izip(cm,rc,pr):
#        print>>output,"#{0}".format("="*20)
#        c.to_csv(output)
#        r.recall.to_csv(output)
#        p.precision.to_csv(output)

def confusion_matrix(contigs,clustering,levels=["family","genus","species"]):
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

if __name__=="__main__":
    parser = ArgumentParser(description="Clustering of metagenomic contigs")
    parser.add_argument('cluster_file', 
        help='Result file from ProBin program (or each line on the form Cluster X, contig1X, contig2X,..)')
    parser.add_argument('phylo_file', 
        help='Phylogenetic file for the cluster_file (each line on the form: contig_id, family, genus, species')
    parser.add_argument('-l', '--level', 
        default='family', type=str, choices=['family','genus','species'],
        help='Calculate statistics for which level.')
    args = parser.parse_args()
    get_statistics(args.cluster_file,args.phylo_file,args.level)