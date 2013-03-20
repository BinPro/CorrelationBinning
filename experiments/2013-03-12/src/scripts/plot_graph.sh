#!/usr/bin/env bash
#SBATCH -A b2010008
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 1-00:00:00

#DATA_PATH="/bubo/home/h20/brynjar/glob/masterproject/DATA/2013-03-05"
DATA_PATH=$HOME"/repos/RESULTS/to_big_for_git/2013-03-08"

#RESULTS_PATH=$HOME"/glob/results/2013-03-08"
RESULTS_PATH=$HOME"/repos/RESULTS/multinomial_kmer_length"

LEVEL="genome"

echo "plotting graph, precision for "$LEVEL" level, contig length 1000, kmer=3-6"
time ../programs/plot_graph.py $DATA_PATH"/plotting_links/kmer-l_3__contig-l_1000" $DATA_PATH"/plotting_links/kmer-l_4__contig-l_1000" $DATA_PATH"/plotting_links/kmer-l_5__contig-l_1000" $DATA_PATH"/plotting_links/kmer-l_6__contig-l_1000" -o $RESULTS_PATH"/figures/multinomial_scores_kmer_3_to_6_1000_"$LEVEL --levels $LEVEL -x "included_contigs_ratio" -y "precision"

echo "plotting graph, precision for "$LEVEL" level, contig length 1000, kmer=7-11"
time ../programs/plot_graph.py $DATA_PATH"/plotting_links/kmer-l_7__contig-l_1000" $DATA_PATH"/plotting_links/kmer-l_8__contig-l_1000" $DATA_PATH"/plotting_links/kmer-l_9__contig-l_1000" $DATA_PATH"/plotting_links/kmer-l_10__contig-l_1000" -o $RESULTS_PATH"/figures/multinomial_scores_kmer_7_to_11_1000_"$LEVEL --levels $LEVEL -x "included_contigs_ratio" -y "precision"
