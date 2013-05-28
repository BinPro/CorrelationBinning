# This experiment aim is to add taxonomic 
# information to the contigs from mock-data.

RESULTS_PATH=$(HOME)"/repos/RESULTS/to_big_for_git/2013-05-02/low_coverage_sample/even/"

parsed_gen_0_0_mock_complete_fixed3.txt:
	ln -s $(RESULTS_PATH)$@ $@

contigs_with_samples_random_abundance_otu.tsv: 
	ln -s $(RESULTS_PATH)$@ $@

add_taxa_to_contigs_time_series.py:
	ln -s ../../src/programs/$@ $@

contigs_with_taxonomic_info.tsv:contigs_with_samples_random_abundance_otu.tsv parsed_gen_0_0_mock_complete_fixed3.txt add_taxa_to_contigs_time_series.py
	./add_taxa_to_contigs_time_series.py $< -o $@ -t parsed_gen_0_0_mock_complete_fixed3.txt
