# Plots all graphs for the master thesis.
RESULTS_PATH="/proj/b2010008/ProBin/RESULTS/"
MULTINOMIAL_PATH=$(RESULTS_PATH)"multinomial_kmer_length/precision_matrix/2013-04-18"


plot_graph.py:
	ln -s '../../../2013-03-12/src/programs/'$@ $@

plot_matrix.py:
	ln -s '../../../2013-03-12/src/programs/'$@ $@

add_rows_to_matrix.sh:
	ln -s '../../src/scripts/'$@ $@

precision_matrix_family.tsv: add_rows_to_matrix.sh 
	./add_rows_to_matrix.sh $(MULTINOMIAL_PATH) family precision_matrix 3 10

precision_matrix_genus.tsv: add_rows_to_matrix.sh
	./add_rows_to_matrix.sh $(MULTINOMIAL_PATH) genus precision_matrix 3 10

precision_matrix_genome.tsv: add_rows_to_matrix.sh
	./add_rows_to_matrix.sh $(MULTINOMIAL_PATH) genome precision_matrix 3 10

multinomial_score_all_kmer_100_family.png: precision_matrix_genome.tsv plot_matrix.py
	./plot_matrix.py $< -o multinomial_score_all_kmer_100 --title 'None'

all: precision_matrix_genome.tsv precision_matrix_genus.tsv precision_matrix_family.tsv

