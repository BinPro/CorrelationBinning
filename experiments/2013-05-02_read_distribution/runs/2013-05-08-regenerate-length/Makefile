
## ##
## Author: Johannes Alneberg 
## e-mail: johannesalneberg@gmail.com
## This makefile executes the experiment to investigate 
## how the names of the strains are used in the fasta file
## and the read mapping file.

FASTA_FILE="/proj/b2010008/ProBin/DATA/chris-mock/references/references.fa"

MASMVALI_FILE="/proj/b2010008/ProBin/DATA/chris-mock/unbalanced/masmvali/masmvaliout/contig-purity.tsv"

PHYLOGENY_FILE="/proj/b2010008/ProBin/DATA/chris-mock/references/phylogeny-references.tsv"

MANUALLY_FIXED_FILE="/proj/b2010008/ProBin/DATA/chris-mock/references/genome_length_fixed.tsv"

genome_length_regenerated.tsv:
	./genome_length.py $(FASTA_FILE) -o genome_length_regenerated.tsv

translated_length.tsv: genome_length_regenerated.tsv
	./translate_column.py genome_length_regenerated.tsv -p $(PHYLOGENY_FILE) --dictwise -o translated_length.tsv

diff_file: translated_length.tsv
	diff $(MANUALLY_FIXED_FILE) translated_length.tsv > diff_file


