#!/usr/bin/env bash
#SBATCH -A b2010008
#SBATCH -p node
#SBATCH -n 8
#SBATCH -t 2-00:00:00
#SBATCH -J Clustering_em

# This experiment aim is to compare how different kmer-lengths affect
# the genomic profile specificity for different taxonomic levels.

DATA_PATH="/proj/b2010008/ProBin/DATA/2013-04-26_clustering_uppmax"
mkdir -p $DATA_PATH

RESULTS_PATH="/proj/b2010008/ProBin/RESULTS/clustering/2013-04-26_clustering_uppmax/compare_kmer_lengths"
mkdir -p $RESULTS_PATH

KMER_LENGTH=" 4 5"
CLUSTER_COUNT="13 55 184"
CLUSTER_ALG="em"

DATA_FILES="contigs_2_2_100_100_start_position.fna contigs_2_2_100_1000_start_position.fna contigs_2_2_100_10000_start_position.fna"
DATE_OUT=$(date +%F_%H:%M)
OUTPUT_ERROR=$RESULTS_PATH/$CLUSTER_ALG"_"$KMER_LENGTH"_"$CLUSTER_COUNT".error_log" 

for ALG in $CLUSTER_ALG
do
  for COUNT in $CLUSTER_COUNT
  do
    for KMER in $KMER_LENGTH
    do
      for INFILE in $DATA_FILES
      do
	OUTPUT_CSV=$RESULTS_PATH/$DATE_OUT"_clustering_"$ALG"_kmer_"$KMER"_clusters_"$COUNT".csv"
	ProBin.py bin $DATA_PATH/$INFILE -k $KMER -c $COUNT -a $ALG -o $OUTPUT_CSV
	wait
      done
    done
  done
done
