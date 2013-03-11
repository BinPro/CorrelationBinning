#!/usr/bin/env bash
#SBATCH -A b2010008
#SBATCH -p node
#SBATCH -t 1-00:00:00

# This experiment aim is to compare how different kmer-lengths affect
# the genomic profile specificity for different taxonomic levels.

DATA_PATH="/bubo/home/h20/brynjar/glob/masterproject/DATA/2013-03-05"
#DATA_PATH=$HOME"/repos/DATA"

RESULTS_PATH=$HOME"/glob/results/2013-03-08"

KMER_LENGTH=7

echo "Scoring contigs, length 100, kmer length "$KMER_LENGTH >&2
time ../programs/score_contigs_multinomial.py $DATA_PATH"/contigs/contigs_2_2_new_100_100.fna" -o $RESULTS_PATH"/score_mul_2_2_100_100_$KMER_LENGTH"".tsv" -t $DATA_PATH"/parsed_taxonomy_2_2_2013-03-05_complete.txt" -k $KMER_LENGTH -d $DATA_PATH"/all_fna"

echo "Scoring contigs, length 1000, kmer length "$KMER_LENGTH >&2
time ../programs/score_contigs_multinomial.py $DATA_PATH"/contigs/contigs_2_2_new_100_1000.fna" -o $RESULTS_PATH"/score_mul_2_2_100_1000_$KMER_LENGTH"".tsv" -t $DATA_PATH"/parsed_taxonomy_2_2_2013-03-05_complete.txt" -k $KMER_LENGTH -d $DATA_PATH"/all_fna"

echo "Scoring contigs, length 10000, kmer_length "$KMER_LENGTH >&2
time ../programs/score_contigs_multinomial.py $DATA_PATH"/contigs/contigs_2_2_new_100_10000.fna" -o $RESULTS_PATH"/score_mul_2_2_100_10000_$KMER_LENGTH"".tsv" -t $DATA_PATH"/parsed_taxonomy_2_2_2013-03-05_complete.txt" -k $KMER_LENGTH -d $DATA_PATH"/all_fna"

