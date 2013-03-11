#!/usr/bin/env bash
#SBATCH -A b2010008
#SBATCH -p node
#SBATCH -t 8:00:00

# This experiment aim is to compare how different kmer-lengths affect
# the genomic profile specificity for different taxonomic levels.

DATA_PATH="/bubo/home/h20/brynjar/glob/masterproject/DATA/2013-03-05"
#DATA_PATH=$HOME"/repos/DATA"

RESULTS_PATH=$HOME"/glob/results/2013-03-08"

#echo "Scoring contigs 3 100" >&2
#time ../programs/score_contigs_multinomial.py $DATA_PATH"/contigs/contigs_2_2_new_100_100.fna" -o $RESULTS_PATH"/score_mul_2_2_100_100.tsv" -t $DATA_PATH"/parsed_taxonomy_2_2_2013-03-05_complete.txt" -k 3 -d $DATA_PATH"/all_fna"

echo "Scoring contigs 3 1000" >&2
time ../programs/score_contigs_multinomial.py $DATA_PATH"/contigs/contigs_2_2_new_100_1000.fna" -o $RESULTS_PATH"/score_mul_2_2_100_1000.tsv" -t $DATA_PATH"/parsed_taxonomy_2_2_2013-03-05_complete.txt" -k 3 -d $DATA_PATH"/all_fna"

echo "Scoring contigs 3 10000" >&2
time ../programs/score_contigs_multinomial.py $DATA_PATH"/contigs/contigs_2_2_new_100_10000.fna" -o $RESULTS_PATH"/score_mul_2_2_100_10000.tsv" -t $DATA_PATH"/parsed_taxonomy_2_2_2013-03-05_complete.txt" -k 3 -d $DATA_PATH"/all_fna"

