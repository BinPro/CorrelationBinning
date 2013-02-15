#!/usr/bin/env bash
#SBATCH -A b2010008
#SBATCH -p node
#SBATCH -t 8:00:00

# This experiment aim is to compare how different kmer-lengths affect
# the genomic profile specificity for different taxonomic levels.

DATA_PATH=$HOME"/glob/data"
#DATA_PATH=$HOME"/repos/DATA"

RESULTS_PATH=$HOME"/glob/results/2013-01-18"

echo "Parsing the taxonomy file" >&2
time ./parse_taxonomy_complete.py -s 2 -g 2 $DATA_PATH"/gen_tax.tab" -o $RESULTS_PATH"/parsed_gen_2_2_u"
echo "Calculate multinomial log probabilities for contigs from all genomes" >&2
echo "k=3"
time ./pairwise_probabilities.py -m "multinomial" -k 3 $RESULTS_PATH"/parsed_gen_2_2_u_complete.txt" -o $RESULTS_PATH"/multinomial_3_1000" -d $DATA_PATH"/reference_genomes_ncbi" -c 1000 -p "genomes" --contig_min_length 1000 --contig_max_length 1000

echo "k=4"
time ./pairwise_probabilities.py -m "multinomial" -k 4 $RESULTS_PATH"/parsed_gen_2_2_u_complete.txt" -o $RESULTS_PATH"/multinomial_4_1000" -d $DATA_PATH"/reference_genomes_ncbi" -c 1000 -p "genomes" --contig_min_length 1000 --contig_max_length 1000

echo "k=5"
time ./pairwise_probabilities.py -m "multinomial" -k 5 $RESULTS_PATH"/parsed_gen_2_2_u_complete.txt" -o $RESULTS_PATH"/multinomial_5_1000" -d $DATA_PATH"/reference_genomes_ncbi" -c 1000 -p "genomes" --contig_min_length 1000 --contig_max_length 1000

echo "Stats summary for k=3" >&2
time ./result_statistics.py $RESULTS_PATH"/multinomial_3_1000" -o $RESULTS_PATH"/multinomial_stats_3_1000"
echo "Stats summary for k=4" >&2
time ./result_statistics.py $RESULTS_PATH"/multinomial_4_1000" -o $RESULTS_PATH"/multinomial_stats_4_1000"
echo "Stats summary for k=5" >&2
time ./result_statistics.py $RESULTS_PATH"/multinomial_5_1000" -o $RESULTS_PATH"/multinomial_stats_5_1000"
