#!/usr/bin/env bash
#SBATCH -A b2010008
#SBATCH -p devel
#SBATCH -t 1:00:00
#SBATCH -J Clustering_kmeans

# This experiment aim is to compare how different kmer-lengths affect
# the genomic profile specificity for different taxonomic levels.

DATA_PATH="/proj/b2010008/ProBin/DATA/2013-04-26"
#DATA_PATH=$HOME"/repos/DATA"

RESULTS_PATH="/proj/b2010008/ProBin/RESULTS/clustering/2013-04-26"
mkdir -p $RESULTS_PATH

KMER_LENGTH=4
CLUSTER_COUNT=7
CLUSTER_ALG="kmeans"

OUTPUT_CSV=$RESULTS_PATH/$CLUSTER_ALG"_"$KMER_LENGTH"_"$CLUSTER_COUNT".csv" 
OUTPUT_ERROR=$RESULTS_PATH/$CLUSTER_ALG"_"$KMER_LENGTH"_"$CLUSTER_COUNT".error_log" 

COMMAND="ProBin.py bin $DATA_PATH/generated_contigs_10000_test.fna -k $KMER_LENGTH -c $CLUSTER_COUNT -a $CLUSTER_ALG -o $OUTPUT_CSV"

echo $COMMAND
time $COMMAND

