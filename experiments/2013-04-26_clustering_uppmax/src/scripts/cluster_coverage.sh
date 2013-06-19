#!/usr/bin/env bash
#SBATCH -A b2010008
#SBATCH -p node
##SBATCH -C mem72GB
#SBATCH -n 8
#SBATCH -t 04:00:00
#SBATCH -J Clustering_coverage
#SBATCH --mail-type=ALL
#SBATCH --mail-user=brynjar.bjarnason@scilifelab.se

DATA_PATH="/proj/b2010008/ProBin/RESULTS/clustering/2013-04-26_clustering_uppmax/mock-even-low-coverage"
RESULTS_PATH="/proj/b2010008/ProBin/RESULTS/clustering/2013-04-26_clustering_uppmax/mock-even-low-coverage/"

CLUSTER_COUNT="40"
CLUSTER_ALG="kmeans"

DATA_FILES="contigs_with_taxonomic_info.tsv"
DATE_OUT=$(date +%F_%H:%M)
OUTPUT_ERROR=$RESULTS_PATH/$CLUSTER_ALG"_"$KMER_LENGTH"_"$CLUSTER_COUNT".error_log" 

for ALG in $CLUSTER_ALG
do
  for COUNT in $CLUSTER_COUNT
  do
    for INFILE in $DATA_FILES
    do
	ProBin.py bin coverage $DATA_PATH/$INFILE 2012-03-25 2013-01-18 -c $COUNT -a $ALG -o $RESULTS_PATH
	wait
    done
  done
done
