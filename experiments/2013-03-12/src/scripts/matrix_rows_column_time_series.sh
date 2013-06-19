#!/usr/bin/env bash
#SBATCH -A b2010008
#SBATCH -p devel
#SBATCH -n 1
#SBATCH -t 1:00:00

#DATA_PATH="/bubo/home/h20/brynjar/glob/masterproject/DATA/2013-03-05"
#DATA_PATH=$HOME"/repos/RESULTS/to_big_for_git/2013-03-11/2013-03-26"
DATA_PATH="/proj/b2010008/ProBin/RESULTS/to_big_for_git/2013-05-17/2013-06-02"

#RESULTS_PATH=$HOME"/glob/results/2013-03-08"
RESULTS_PATH="/proj/b2010008/ProBin/RESULTS/mock_even_isotropic/precision_matrix/2013-06-02"

LEVEL=$1
KMER_LENGTH=$2
time ../programs/matrix_row.py $DATA_PATH"/score_time_series_even_only_iso_5_columns.tsv" -o $RESULTS_PATH"/m_row_even_40_2_5_"$LEVEL --level $LEVEL -f "precision"

time ../programs/matrix_row.py $DATA_PATH"/score_time_series_even_only_iso_10_columns.tsv" -o $RESULTS_PATH"/m_row_even_40_2_10_"$LEVEL --level $LEVEL -f "precision"
