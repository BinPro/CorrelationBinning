#!/usr/bin/env bash
#SBATCH -A b2010008
#SBATCH -p devel
#SBATCH -n 1
#SBATCH -t 1:00:00

DATA_PATH="/proj/b2010008/ProBin/RESULTS/to_big_for_git/2013-05-17/2013-05-28"

RESULTS_PATH="/proj/b2010008/ProBin/RESULTS/mock_even_isotropic/precision_matrix/2013-05-28"

LEVEL=$1
KMER_LENGTH=$2
time ../programs/matrix_row.py $DATA_PATH"/score_time_series_even_only_iso_2_columns.tsv" -o $RESULTS_PATH"/m_row_even_40_2_"$KMER_LENGTH"_"$LEVEL --level $LEVEL -f "precision"