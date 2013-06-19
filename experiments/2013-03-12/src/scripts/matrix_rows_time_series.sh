#!/usr/bin/env bash
#SBATCH -A b2010008
#SBATCH -p devel
#SBATCH -n 1
#SBATCH -t 1:00:00

DATA_PATH="/proj/b2010008/ProBin/RESULTS/to_big_for_git/2013-05-17/2013-06-04"

RESULTS_PATH="/proj/b2010008/ProBin/RESULTS/mock_even_isotropic/precision_matrix/2013-06-04"

LEVEL=$1
KMER_LENGTH=$2
time ../programs/matrix_row.py $DATA_PATH"/score_time_series_real_only_iso.tsv" -o $RESULTS_PATH"/m_row_even_40_real_"$KMER_LENGTH"_"$LEVEL --level $LEVEL -f "precision"

time ../programs/matrix_row.py $DATA_PATH"/score_time_series_cv_only_iso.tsv" -o $RESULTS_PATH"/m_row_even_40_cv_"$KMER_LENGTH"_"$LEVEL --level $LEVEL -f "precision"