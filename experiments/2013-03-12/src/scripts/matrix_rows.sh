#!/usr/bin/env bash
#SBATCH -A b2010008
#SBATCH -p core
#SBATCH -t 1:00:00

DATA_PATH="/proj/b2010008/ProBin/RESULTS/to_big_for_git/2013-03-08/2013-04-18"

RESULTS_PATH="/proj/b2010008/ProBin/RESULTS/multinomial_kmer_length/precision_matrix/2013-04-18"

LEVEL=$1
KMER_MIN=$2
KMER_MAX=$3
KMER_LENGTH=$KMER_MIN
while [[ $KMER_LENGTH -le $KMER_MAX ]]
do
    echo "calculating matrix row, precision for "$LEVEL" level, contig length 100,1000 and 10000, kmer="$KMER_LENGTH
    time ./matrix_row.py $DATA_PATH"/score_mul_2_2_100_100_"$KMER_LENGTH".tsv" $DATA_PATH"/score_mul_2_2_100_1000_"$KMER_LENGTH".tsv" $DATA_PATH"/score_mul_2_2_100_10000_"$KMER_LENGTH".tsv"  -o $RESULTS_PATH"/m_row_100_1000_10000_100_"$KMER_LENGTH"_"$LEVEL --level $LEVEL -f "precision"
    ((KMER_LENGTH = KMER_LENGTH + 1))
done