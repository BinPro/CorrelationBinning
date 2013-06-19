#!/usr/bin/env bash
#SBATCH -A b2010008
#SBATCH -p devel
#SBATCH -n 1
#SBATCH -t 1:00:00

DATA_PATH="/proj/b2010008/ProBin/RESULTS/to_big_for_git/2013-06-04/2013-06-04"

RESULTS_PATH="/proj/b2010008/ProBin/RESULTS/simple_add/precision_matrix/2013-06-04"

LEVEL=$1
KMER_LENGTH=$2
time ../programs/matrix_row.py $DATA_PATH"/score_combined_real_46_columns_"$KMER_LENGTH"_100.tsv" -o $RESULTS_PATH"/m_row_100_real_40_"$KMER_LENGTH"_"$LEVEL --level $LEVEL -f "precision"