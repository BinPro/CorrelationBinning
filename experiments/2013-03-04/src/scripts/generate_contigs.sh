#!/usr/bin/env bash
#SBATCH -A b2010008
#SBATCH -p node
#SBATCH -t 4:00:00

# This script generates contigs to be used for classification 
# evaluation. The output is saved in order to be reused for all 
# experiments.

#DATA_PATH=$HOME"/glob/data"
DATA_PATH=$HOME"/repos/DATA"

#RESULTS_PATH=$HOME"/glob/results/2013-01-18"
RESULTS_PATH=$HOME"/repos/DATA/contigs"

echo "Generate contigs" >&2
time ../programs/generate_contigs.py $DATA_PATH"/parsed_gen_2_2_complete_old.txt" -o $RESULTS_PATH"/contigs_2_2_old_100_100.fna" -d $DATA_PATH"/reference_genomes_ncbi" -c 100 -p "genomes" --contig_min_length 100 --contig_max_length 100
