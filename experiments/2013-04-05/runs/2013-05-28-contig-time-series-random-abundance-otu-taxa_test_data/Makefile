# This experiment aim is to add taxonomic 
# information to the contigs from mock-data.

DATA_PATH=$(HOME)"/repos/Correlation-Binning/src/tests/fixtures/"


parsed_gen_0_0_mock_test_complete_2013-05-24.txt:
	ln -s $(DATA_PATH)$@ $@

contigs_with_samples_otu_test.tsv: 
	ln -s $(DATA_PATH)$@ $@

add_taxa_to_contigs_time_series.py:
	ln -s ../../src/programs/$@ $@

contigs_with_taxonomic_info.tsv:contigs_with_samples_otu_test.tsv parsed_gen_0_0_mock_test_complete_2013-05-24.txt add_taxa_to_contigs_time_series.py
	./add_taxa_to_contigs_time_series.py $< -o $@ -t parsed_gen_0_0_mock_test_complete_2013-05-24.txt
