#!/usr/bin/env bash
#SBATCH -A b2010008
#SBATCH -p node
#SBATCH -t 3-00:00:00

# This experiment aim is to compare how different kmer-lengths affect
# the genomic profile specificity for different taxonomic levels.

DATA_PATH="/bubo/home/h20/brynjar/glob/masterproject/DATA/2013-03-05"
#DATA_PATH=$HOME"/repos/DATA"
#DATA_PATH=$HOME"/glob/RESULTS"

RESULTS_PATH=$HOME"/glob/RESULTS/to_big_for_git/2013-03-11/2013-04-10"
#RESULTS_PATH=$HOME"/glob/results/2013-03-08/2013-03-18"

KMER_LENGTH=10

CONTIG_LENGTH=100
echo "Scoring contigs, length "$CONTIG_LENGTH", kmer_length "$KMER_LENGTH >&2
time ../programs/score_contigs_dirichlet.py $DATA_PATH"/contigs/start_position/contigs_2_2_100_"$CONTIG_LENGTH"_start_position.fna" -o $RESULTS_PATH"/score_dirichlet_2_2_100_"$CONTIG_LENGTH"_"$KMER_LENGTH".tsv" -t $DATA_PATH"/parsed_taxonomy_2_2_2013-03-05_complete.txt" -k $KMER_LENGTH -d $DATA_PATH"/all_fna" -c $CONTIG_LENGTH

CONTIG_LENGTH=1000
echo "Scoring contigs, length "$CONTIG_LENGTH", kmer_length "$KMER_LENGTH >&2
time ../programs/score_contigs_dirichlet.py $DATA_PATH"/contigs/start_position/contigs_2_2_100_"$CONTIG_LENGTH"_start_position.fna" -o $RESULTS_PATH"/score_dirichlet_2_2_100_"$CONTIG_LENGTH"_"$KMER_LENGTH".tsv" -t $DATA_PATH"/parsed_taxonomy_2_2_2013-03-05_complete.txt" -k $KMER_LENGTH -d $DATA_PATH"/all_fna" -c $CONTIG_LENGTH

CONTIG_LENGTH=10000
echo "Scoring contigs, length "$CONTIG_LENGTH", kmer_length "$KMER_LENGTH >&2
time ../programs/score_contigs_dirichlet.py $DATA_PATH"/contigs/start_position/contigs_2_2_100_"$CONTIG_LENGTH"_start_position.fna" -o $RESULTS_PATH"/score_dirichlet_2_2_100_"$CONTIG_LENGTH"_"$KMER_LENGTH".tsv" -t $DATA_PATH"/parsed_taxonomy_2_2_2013-03-05_complete.txt" -k $KMER_LENGTH -d $DATA_PATH"/all_fna" -c $CONTIG_LENGTH

