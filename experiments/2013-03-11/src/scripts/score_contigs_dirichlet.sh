#!/usr/bin/env bash
#SBATCH -A b2010008
#SBATCH -p devel
#SBATCH -t 59:00

# This experiment aim is to compare how different kmer-lengths affect
# the genomic profile specificity for different taxonomic levels.

DATA_PATH="/bubo/home/h20/brynjar/glob/masterproject/DATA/2013-03-05"
#DATA_PATH=$HOME"/repos/DATA"
#DATA_PATH=$HOME"/glob/RESULTS"

RESULTS_PATH=$HOME"/glob/RESULTS/to_big_for_git/2013-03-11/2013-03-26"
#RESULTS_PATH=$HOME"/glob/results/2013-03-08/2013-03-18"

KMER_LENGTH=3

#echo "Scoring contigs, length 100, kmer length "$KMER_LENGTH >&2
#time ../programs/score_contigs_dirichlet.py $DATA_PATH"/contigs/start_position/contigs_2_2_100_100_start_position.fna" -o $RESULTS_PATH"/score_dirichlet_2_2_100_100_$KMER_LENGTH"".tsv" -t $DATA_PATH"/parsed_taxonomy_2_2_2013-03-05_complete.txt" -k $KMER_LENGTH -d $DATA_PATH"/all_fna"

#echo "Scoring contigs, length 1000, kmer length "$KMER_LENGTH >&2
#time ../programs/score_contigs_dirichlet.py $DATA_PATH"/contigs/start_position/contigs_2_2_100_1000_start_position.fna" -o $RESULTS_PATH"/score_dirichlet_2_2_100_1000_$KMER_LENGTH"".tsv" -t $DATA_PATH"/parsed_taxonomy_2_2_2013-03-05_complete.txt" -k $KMER_LENGTH -d $DATA_PATH"/all_fna"

CONTIG_LENGTH=10000
echo "Scoring contigs, length "$CONTIG_LENGTH", kmer_length "$KMER_LENGTH >&2
time ../programs/score_contigs_dirichlet.py $DATA_PATH"/contigs/start_position/contigs_2_2_100_"$CONTIG_LENGTH"_start_pos_test.fna" -o $RESULTS_PATH"/score_dirichlet_2_2_100_"$CONTIG_LENGTH"_"$KMER_LENGTH"_test.tsv" -t $DATA_PATH"/contigs/parsed_gen_2_2_test.txt" -k $KMER_LENGTH -d $DATA_PATH"/all_fna" -c $CONTIG_LENGTH

