#!/usr/bin/env bash
#SBATCH -A b2010008
#SBATCH -p node
#SBATCH -t 8:00:00

# This experiment aim is to compare how different kmer-lengths affect
# the genomic profile specificity for different taxonomic levels.

#DATA_PATH=$HOME"/glob/data"
DATA_PATH=$HOME"/repos/DATA"

RESULTS_PATH=$HOME"/glob/results/2013-03-08"

echo "Scoring contigs" >&2
time ../programs/score_contigs_multinomial.py $DATA_PATH"/parsed