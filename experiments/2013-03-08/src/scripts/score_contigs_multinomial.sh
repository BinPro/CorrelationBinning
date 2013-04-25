#!/usr/bin/env bash
#SBATCH -A b2010008
#SBATCH -p devel
#SBATCH -t 1:00:00

# This experiment aim is to compare how different kmer-lengths affect
# the genomic profile specificity for different taxonomic levels.

DATA_PATH="/proj/b2010008/ProBin/DATA/2013-03-05"
#DATA_PATH=$HOME"/repos/DATA"

RESULTS_PATH="/proj/b2010008/ProBin/RESULTS/to_big_for_git/2013-03-08/2013-04-18"

KMER_LENGTH=2

for CONTIG_LENGTH in 100 1000
do
    echo "Scoring contigs, length "$CONTIG_LENGTH", kmer length "$KMER_LENGTH >&2
    time ../programs/score_contigs_multinomial.py $DATA_PATH"/contigs/start_position/contigs_2_2_100_"$CONTIG_LENGTH"_start_position.fna" -o $RESULTS_PATH"/score_mul_2_2_100_"$CONTIG_LENGTH"_"$KMER_LENGTH".tsv" -t $DATA_PATH"/parsed_taxonomy_2_2_2013-03-05_complete.txt" -k $KMER_LENGTH -d $DATA_PATH"/all_fna"
done