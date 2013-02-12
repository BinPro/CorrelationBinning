#!/usr/bin/env bash
# Run me like this: ./test_bash.sh > stdout.txt 2> stderr.txt
# This experiment aim is to compare how different kmer-lengths affect 
# the genomic profile specificity for different taxonomic levels.
echo "Parsing the taxonomy file" >&2
#time ./parse_taxonomy_complete.py -s 2 -g 2 ../data/gen_tax.tab -o ../results/parsed_gen_2_2
echo "Calculate multinomial log probabilities for contigs from all genomes" >&2
time ./pairwise_probabilities.py -m "multinomial" -k 4 ../results/parsed_gen_2_2_test.txt -o ../results/multinomial_prob_4_test.txt -d /home/johannes/repos/DATA/reference_genomes_ncbi -c 100 -p "genomes" --contig_min_length 1000 --contig_max_length 1000

echo "Stats summary for k=4" >&2
time ./result_statistics.py ../results/multinomial_prob_families.txt -o ../results/multinomial_stats_4.txt
