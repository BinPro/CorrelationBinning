## ##
## Author: Johannes Alneberg 
## e-mail: johannesalneberg@gmail.com
## This makefile executes the experiment of read mapping
## to compare the distribution of reads for contigs from 
## the same genome.

CONTIGS_FILE="/proj/b2010008/ProBin/DATA/chris-mock/velvet-noscaf-newbler/even-contigs.fa"

MASMVALI_FILE="/proj/b2010008/ProBin/DATA/chris-mock/even/masmvali/masmvaliout/contig-purity.tsv"



full_read_mappings.csv:
	ln -s $(MASMVALI_FILE) full_read_mappings.csv

contigs.fna:
	ln -s $(CONTIGS_FILE) contigs.fna

read_mappings.csv: full_read_mappings.csv contigs.fna
	ProBin.py preprocess full_read_mappings.csv -f masmvali -o read_mappings.csv -c contigs.fna -s

