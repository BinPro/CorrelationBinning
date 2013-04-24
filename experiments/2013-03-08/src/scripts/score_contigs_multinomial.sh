#!/usr/bin/env bash
#SBATCH -A b2010008
#SBATCH -p node
#SBATCH -t 1-00:00:00

# This experiment aim is to compare how different kmer-lengths affect
# the genomic profile specificity for different taxonomic levels.

DATA_PATH="/glob/brynjar/masterproject/DATA/2013-03-05"
#DATA_PATH=$HOME"/repos/DATA"

RESULTS_PATH=$HOME"/glob/RESULTS/to_big_for_git/2013-03-08/2013-04-18"

KMER_LENGTH=11

CONTIG_LENGTH=100
echo "Scoring contigs, length "$CONTIG_LENGTH", kmer length "$KMER_LENGTH >&2
time ../programs/score_contigs_multinomial.py $DATA_PATH"/contigs/start_position/contigs_2_2_100_"$CONTIG_LENGTH"_start_position.fna" -o $RESULTS_PATH"/score_mul_2_2_100_"$CONTIG_LENGTH"_"$KMER_LENGTH".tsv" -t $DATA_PATH"/parsed_taxonomy_2_2_2013-03-05_complete.txt" -k $KMER_LENGTH -d $DATA_PATH"/all_fna"

CONTIG_LENGTH=1000
echo "Scoring contigs, length 1000, kmer length "$KMER_LENGTH >&2
time ../programs/score_contigs_multinomial.py $DATA_PATH"/contigs/start_position/contigs_2_2_100_"$CONTIG_LENGTH"_start_position.fna" -o $RESULTS_PATH"/score_mul_2_2_100_"$CONTIG_LENGTH"_"$KMER_LENGTH".tsv" -t $DATA_PATH"/parsed_taxonomy_2_2_2013-03-05_complete.txt" -k $KMER_LENGTH -d $DATA_PATH"/all_fna"

CONTIG_LENGTH=10000
echo "Scoring contigs, length 10000, kmer_length "$KMER_LENGTH >&2
time ../programs/score_contigs_multinomial.py $DATA_PATH"/contigs/start_position/contigs_2_2_100_"$CONTIG_LENGTH"_start_position.fna" -o $RESULTS_PATH"/score_mul_2_2_100_"$CONTIG_LENGTH"_"$KMER_LENGTH".tsv" -t $DATA_PATH"/parsed_taxonomy_2_2_2013-03-05_complete.txt" -k $KMER_LENGTH -d $DATA_PATH"/all_fna"

