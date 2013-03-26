#!/usr/bin/env bash
#SBATCH -A b2010008
#SBATCH -p node
#SBATCH -t 4:00:00

# This script generates genome parts to be used for classification 
# evaluation. The output is saved in order to be reused for all 
# experiments.

DATA_PATH="/bubo/home/h20/brynjar/glob/masterproject/DATA/2013-03-05"

#DATA_PATH=$HOME"/repos/DATA"

RESULTS_PATH=$HOME"/glob/RESULTS/to_big_for_git/genome_parts/"
#RESULTS_PATH=$HOME"/repos/DATA/contigs"

echo "Generate genome parts" >&2
time ../programs/split_genomes.py $DATA_PATH"/parsed_taxonomy_2_2_2013-03-05_complete.txt" -o $RESULTS_PATH"/genome_parts_2_2_new_100_10000.fna" -d $DATA_PATH"/all_fna" --part_length 10000
