#!/usr/bin/env bash
#SBATCH -A b2010008
#SBATCH -p node
#SBATCH -t 4:00:00

# This script generates contigs to be used for classification 
# evaluation. The output is saved in order to be reused for all 
# experiments.

#DATA_PATH="/bubo/home/h20/brynjar/glob/masterproject/DATA/2013-03-05"
DATA_PATH=$HOME"/repos/DATA/2013-03-05"

#RESULTS_PATH="/bubo/home/h20/brynjar/glob/masterproject/DATA/2013-03-05/contigs"
RESULTS_PATH=$HOME"/repos/RESULTS/contigs/start_position"

echo "Generate contigs" >&2
time ../programs/generate_contigs.py $DATA_PATH"/parsed_taxonomy_2_2_2013-03-05_complete.txt" -o $RESULTS_PATH"/contigs_2_2_100_100_start_position.fna" -d $DATA_PATH"/all_fna" -c 100 -p "genomes" --contig_min_length 100 --contig_max_length 100 --start_position
