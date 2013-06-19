#!/usr/bin/env bash
#SBATCH -A b2010008
#SBATCH -p devel
#SBATCH -n 1
#SBATCH -t 1:00:00

#DATA_PATH="/bubo/home/h20/brynjar/glob/masterproject/DATA/2013-03-05"
#DATA_PATH=$HOME"/repos/RESULTS/to_big_for_git/2013-03-11/2013-03-26"
DATA_PATH="/proj/b2010008/ProBin/RESULTS/to_big_for_git/2013-03-08/2013-06-05"

#RESULTS_PATH=$HOME"/glob/results/2013-03-08"
RESULTS_PATH="/proj/b2010008/ProBin/RESULTS/multinomial_kmer_length/precision_matrix/2013-06-03"

LEVEL=$1
KMER_MIN=$2
KMER_MAX=$3
KMER_LENGTH=$KMER_MIN
while [[ $KMER_LENGTH -le $KMER_MAX ]]
do
    echo "calculating matrix row, precision for "$LEVEL" level, contig length 100,1000 and 10000, kmer="$KMER_LENGTH
    time ../programs/matrix_row.py $DATA_PATH"/score_mul_0_0_real_"$KMER_LENGTH".tsv" -o $RESULTS_PATH"/m_row_real_"$KMER_LENGTH"_"$LEVEL --level $LEVEL -f "precision"
    ((KMER_LENGTH = KMER_LENGTH + 1))
done