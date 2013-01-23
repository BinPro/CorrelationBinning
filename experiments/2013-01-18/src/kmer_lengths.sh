#!/usr/bin/env bash
# This experiment aim is to compare how different kmer-lengths affect 
# the genomic profile specificity for different taxonomic levels.
echo "Parsing the taxonomy file"
./parse_taxonomy.py -s 4 -g 4 ../data/gen_tax.tab -o ../results/parsed_gen
echo "Comparing families"
./pairwise_probabilities.py -m "multinomial" -k 4 ../results/parsed_gen_families.txt -o ../results/multinomial_prob_families.txt
echo "Comparing genera"
./pairwise_probabilities.py -m "multinomial" -k 4 ../results/parsed_gen_genera.txt -o ../results/multinomial_prob_genera.txt
