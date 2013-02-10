#!/usr/bin/env bash
# This experiment aim is to compare how different kmer-lengths affect 
# the genomic profile specificity for different taxonomic levels.
echo "Parsing the taxonomy file"
./parse_taxonomy_complete.py -s 2 -g 2 ../data/gen_tax.tab -o ../results/parsed_gen_2_2
echo "Calculate multinomial log probabilities for contigs from all genomes"
./pairwise_probabilities.py -m "multinomial" -k 4 ../results/parsed_gen_2_2_complete.txt -o ../results/multinomial_prob_families.txt -d /home/johannes/repos/DATA/reference_genomes_ncbi -c 100 -p "genomes" --contig_min_length 1000 --contig_max_length 1000
