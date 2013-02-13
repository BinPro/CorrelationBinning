#!/usr/bin/env bash
#SBATCH -A b2010008
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 20:00

# Run me like this: ./test_bash.sh > stdout.txt 2> stderr.txt
# This experiment aim is to compare how different kmer-lengths affect

DATA_PATH=$HOME"/glob/data"
#DATA_PATH=$HOME"/repos/DATA/reference_genomes_ncbi"

# the genomic profile specificity for different taxonomic levels.
echo "Parsing the taxonomy file" >&2
time ./parse_taxonomy_complete.py -s 2 -g 2 ../data/gen_tax.tab -o ../results/parsed_gen_2_2

echo "Calculate multinomial log probabilities for contigs from all genomes" >&2
echo "k=3"
time ./pairwise_probabilities.py -m "multinomial" -k 3 ../results/parsed_gen_2_2_test.txt -o ../results/upp_test_3 -d $DATA_PATH -c 100 -p "genomes" --contig_min_length 1000 --contig_max_length 1000

echo "k=4" >&2
time ./pairwise_probabilities.py -m "multinomial" -k 4 ../results/parsed_gen_2_2_test.txt -o ../results/upp_test_4 -d $DATA_PATH -c 100 -p "genomes" --contig_min_length 1000 --contig_max_length 1000

echo "k=5" >&2
time ./pairwise_probabilities.py -m "multinomial" -k 5 ../results/parsed_gen_2_2_test.txt -o ../results/upp_test_5 -d $DATA_PATH -c 100 -p "genomes" --contig_min_length 1000 --contig_max_length 1000

echo "Stats summary for k=3" >&2
time ./result_statistics.py ../results/upp_test_3 -o ../results/upp_test_stat_3

echo "Stats summary for k=4" >&2
time ./result_statistics.py ../results/upp_test_4 -o ../results/upp_test_stat_4

echo "Stats summary for k=5" >&2
time ./result_statistics.py ../results/upp_test_5 -o ../results/upp_test_stat_5
