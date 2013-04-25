#!/usr/bin/env bash
#SBATCH -A b2010008
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 1-00:00:00

#DATA_PATH="/bubo/home/h20/brynjar/glob/masterproject/DATA/2013-03-05"
#DATA_PATH=$HOME"/repos/RESULTS/to_big_for_git/2013-03-11/2013-03-26"
DATA_PATH=$HOME"/glob/RESULTS/to_big_for_git/2013-03-11/2013-04-18"

#RESULTS_PATH=$HOME"/glob/results/2013-03-08"
RESULTS_PATH=$HOME"/glob/RESULTS/dirichlet_multinomial_kmer_length/precision_matrix/2013-04-18"

LEVEL=$1
KMER_MIN=$2
KMER_MAX=$3
KMER_LENGTH=$KMER_MIN
while [[ $KMER_LENGTH -le $KMER_MAX ]]
do
    echo "calculating matrix row, precision for "$LEVEL" level, contig length 100,1000 and 10000, kmer="$KMER_LENGTH
    time ../programs/matrix_row.py $DATA_PATH"/score_dirichlet_2_2_100_100_"$KMER_LENGTH".tsv" $DATA_PATH"/score_dirichlet_2_2_100_1000_"$KMER_LENGTH".tsv" $DATA_PATH"/score_dirichlet_2_2_100_10000_"$KMER_LENGTH".tsv"  -o $RESULTS_PATH"/m_row_100_1000_10000_100_"$KMER_LENGTH"_"$LEVEL --level $LEVEL -f "precision"
    ((KMER_LENGTH = KMER_LENGTH + 1))
done