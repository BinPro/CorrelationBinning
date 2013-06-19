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

precision_matrix_species.tsv: add_rows_to_matrix.sh
	./add_rows_to_matrix.sh $(MULTINOMIAL_PATH) species precision_matrix 3 10

multinomial_score_all_kmer_100_family.png: precision_matrix_family.tsv plot_matrix.py
	./plot_matrix.py $< -o multinomial_score_all_kmer_100_family --title 'Multinomial Family Precision'

multinomial_score_all_kmer_100_genus.png: precision_matrix_genus.tsv plot_matrix.py
	./plot_matrix.py $< -o multinomial_score_all_kmer_100_genus --title 'Multinomial Genus Precision'

multinomial_score_all_kmer_100_species.png: precision_matrix_species.tsv plot_matrix.py
	./plot_matrix.py $< -o multinomial_score_all_kmer_100_species --title 'Multinomial Species Precision'


all:multinomial_score_all_kmer_100_species.png multinomial_score_all_kmer_100_genus.png multinomial_score_all_kmer_100_family.png
	echo "Done, man!"
	touch all
