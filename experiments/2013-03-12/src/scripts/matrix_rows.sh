#!/usr/bin/env bash
#SBATCH -A b2010008
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 1-00:00:00

#DATA_PATH="/bubo/home/h20/brynjar/glob/masterproject/DATA/2013-03-05"
DATA_PATH=$HOME"/repos/RESULTS/to_big_for_git/2013-03-08/links_to_best"

#RESULTS_PATH=$HOME"/glob/results/2013-03-08"
RESULTS_PATH=$HOME"/repos/RESULTS/multinomial_kmer_length/precision_matrix/2013-03-20"

LEVEL=$1
KMER_LENGTH=$2

echo "calculating matrix row, precision for "$LEVEL" level, contig length 100,1000 and 10000, kmer="$KMER_LENGTH
time ../programs/matrix_row.py $DATA_PATH"/score_mul_2_2_100_100_"$KMER_LENGTH".tsv" $DATA_PATH"/score_mul_2_2_100_1000_"$KMER_LENGTH".tsv" $DATA_PATH"/score_mul_2_2_100_10000_"$KMER_LENGTH".tsv"  -o $RESULTS_PATH"/m_row_100_1000_10000_100_"$KMER_LENGTH"_"$LEVEL --level $LEVEL -f "precision"
