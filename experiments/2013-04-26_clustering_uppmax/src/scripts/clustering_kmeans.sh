#!/usr/bin/env bash
#SBATCH -A b2010008
#SBATCH -p node
#SBATCH -t 1-00:00:00
#SBATCH -J Clustering_kmeans

# This experiment aim is to compare how different kmer-lengths affect
# the genomic profile specificity for different taxonomic levels.

DATA_PATH="/proj/b2010008/ProBin/DATA/2013-04-26"
#DATA_PATH=$HOME"/repos/DATA"

RESULTS_PATH=$HOME"/proj/b2010008/ProBin/RESULTS/2013-04-26"
mkdir -p $RESULTS_PATH

KMER_LENGTH=3
CLUSTER_COUNT=3
CLUSTER_ALG="kmeans"

COMMAND=ProBin.py $DATA_PATH/........ -k $KMER_LENGTH -c $CLUSTER_COUNT -a $CLUSTER_ALG -o $RESULT_PATH/$CLUSTER_ALG"_"$KMER_LENGTH"_"$CLUSTER_COUNT".csv" 

echo $COMMAND
time $COMMAND

