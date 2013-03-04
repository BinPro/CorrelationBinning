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
RESULTS_PATH=$HOME"/repos/Correlation-Binning/experiments/2013-03-04/results"

echo "Parsing the taxonomy file" >&2
time ../programs/parse_taxonomy_complete.py -s 2 -g 2 $DATA_PATH"/gen_tax.tab" -o $RESULTS_PATH"/parsed_gen_2_2_u"

echo "Generate contigs" >&2
time ../programs/generate_contigs.py $RESULTS_PATH"/parsed_gen_2_2_u" -o $RESULTS_PATH"/multinomial_3_1000" -d $DATA_PATH"/reference_genomes_ncbi" -c 1000 -p "genomes" --contig_min_length 1000 --contig_max_length 1000
