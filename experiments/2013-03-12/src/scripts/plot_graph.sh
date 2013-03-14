#!/usr/bin/env bash
#SBATCH -A b2010008
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 1-00:00:00

DATA_PATH="/bubo/home/h20/brynjar/glob/masterproject/DATA/2013-03-05"
#DATA_PATH=$HOME"/repos/DATA"

RESULTS_PATH=$HOME"/glob/results/2013-03-08"
echo "plotting graph, precision for genus level, contig length 1000"
time ../programs/plot_graph.py $RESULTS_PATH"/plotting_links/kmer-l_3_contig-l_1000" $RESULTS_PATH"/plotting_links/kmer-l_4_contig-l_1000" $RESULTS_PATH"/plotting_links/kmer-l_5_contig-l_1000" $RESULTS_PATH"/plotting_links/kmer-l_6_contig-l_1000" $RESULTS_PATH"/plotting_links/kmer-l_7_contig-l_1000" $RESULTS_PATH"/plotting_links/kmer-l_8_contig-l_1000" $RESULTS_PATH"/plotting_links/kmer-l_9_contig-l_1000" $RESULTS_PATH"/plotting_links/kmer-l_10_contig-l_1000" -o $RESULTS_PATH"/figures/multinomial_scores_all_kmer_1000_genus" --levels "genus" -x "included_contigs_ratio" -y "precision"

echo "plotting graph, precision for genome level, contig length 1000"
time ../programs/plot_graph.py $RESULTS_PATH"/plotting_links/kmer-l_3_contig-l_1000" $RESULTS_PATH"/plotting_links/kmer-l_4_contig-l_1000" $RESULTS_PATH"/plotting_links/kmer-l_5_contig-l_1000" $RESULTS_PATH"/plotting_links/kmer-l_6_contig-l_1000" $RESULTS_PATH"/plotting_links/kmer-l_7_contig-l_1000" $RESULTS_PATH"/plotting_links/kmer-l_8_contig-l_1000" $RESULTS_PATH"/plotting_links/kmer-l_9_contig-l_1000" $RESULTS_PATH"/plotting_links/kmer-l_10_contig-l_1000" -o $RESULTS_PATH"/figures/multinomial_scores_all_kmer_1000_genome" --levels "genome" -x "included_contigs_ratio" -y "precision"
