# Plots all graphs for the master thesis.
RESULTS_PATH=$(HOME)"/repos/RESULTS/"
MULTINOMIAL_PATH=$(RESULTS_PATH)"multinomial_kmer_length/precision_matrix/2013-04-18"
ISOTROPIC_PATH=$(RESULTS_PATH)"mock_even_isotropic/precision_matrix/2013-05-28"


plot_graph.py:
	ln -s '../../../2013-03-12/src/programs/'$@ $@

plot_matrix.py:
	ln -s '../../../2013-03-12/src/programs/'$@ $@

add_rows_to_matrix.sh:
	ln -s '../../src/scripts/'$@ $@

add_rows_to_matrix_mock.sh:
	ln -s '../../src/scripts/'$@ $@

precision_matrix_family.tsv: add_rows_to_matrix.sh 
	./add_rows_to_matrix.sh $(MULTINOMIAL_PATH) family precision_matrix 3 10

precision_matrix_genus.tsv: add_rows_to_matrix.sh
	./add_rows_to_matrix.sh $(MULTINOMIAL_PATH) genus precision_matrix 3 10

precision_matrix_species.tsv: add_rows_to_matrix.sh
	./add_rows_to_matrix.sh $(MULTINOMIAL_PATH) species precision_matrix 3 10

multinomial_in_silico_kmer_length_together.png: precision_matrix_species.tsv precision_matrix_genus.tsv precision_matrix_family.tsv plot_matrix.py
	./plot_matrix.py precision_matrix_species.tsv precision_matrix_genus.tsv precision_matrix_family.tsv -o $@ --title 'Multinomial Together Precision' --together

mock_precision_matrix.tsv: add_rows_to_matrix_mock.sh
	./add_rows_to_matrix_mock.sh $(ISOTROPIC_PATH) species mock_precision_matrix

isotropic_mock_multiple_samples.png: mock_precision_matrix.tsv plot_matrix.py
	./plot_matrix.py $< -o $@ --title 'Mock Multiple Sample Precision' --mock_multiple
all: multinomial_in_silico_kmer_length_together.png isotropic_mock_multiple_samples.png
	echo "Done, man!"
	touch all
