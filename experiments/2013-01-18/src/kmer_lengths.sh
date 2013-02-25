#!/usr/bin/env bash
#SBATCH -A b2010008
#SBATCH -p node
#SBATCH -t 2-00:00:00

# This experiment aim is to compare how different kmer-lengths affect
# the genomic profile specificity for different taxonomic levels.

DATA_PATH=$HOME"/glob/data"
#DATA_PATH=$HOME"/repos/DATA"

RESULTS_PATH=$HOME"/glob/results/2013-01-18"

echo "Parsing the taxonomy file" >&2
time ./parse_taxonomy_complete.py -s 2 -g 2 $DATA_PATH"/gen_tax.tab" -o $RESULTS_PATH"/parsed_gen_2_2_u"
echo "Calculate multinomial log probabilities for contigs from all genomes" >&2
echo "k=6"
time ./pairwise_probabilities.py -m "multinomial" -k 6 $RESULTS_PATH"/parsed_gen_2_2_u_complete.txt" -o $RESULTS_PATH"/multinomial_6_100" -d $DATA_PATH"/reference_genomes_ncbi" -c 100 -p "genomes" --contig_min_length 1000 --contig_max_length 1000

echo "k=7"
time ./pairwise_probabilities.py -m "multinomial" -k 7 $RESULTS_PATH"/parsed_gen_2_2_u_complete.txt" -o $RESULTS_PATH"/multinomial_7_100" -d $DATA_PATH"/reference_genomes_ncbi" -c 100 -p "genomes" --contig_min_length 1000 --contig_max_length 1000

echo "Stats summary for k=6" >&2
time ./result_statistics.py $RESULTS_PATH"/multinomial_6_100" -o $RESULTS_PATH"/multinomial_stats_6_100"
echo "Stats summary for k=7" >&2
time ./result_statistics.py $RESULTS_PATH"/multinomial_7_100" -o $RESULTS_PATH"/multinomial_stats_7_100"
