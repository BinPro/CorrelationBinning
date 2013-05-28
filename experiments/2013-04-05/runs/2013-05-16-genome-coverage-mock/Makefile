# This experiment aim is to trim the kalmar time-series file
# to the number of genomes present in the mock community

RESULTS_PATH=$(HOME)"/repos/RESULTS/to_big_for_git/2013-05-02/low_coverage_sample/unbalanced/"

calculate_genome_coverage.py:
	ln -s "../../src/programs/"$@ $@

read_mappings_with_coverage.tsv: calculate_genome_coverage.py
	time ./$< $(RESULTS_PATH)"read_mappings.csv" --read_mapping_column "full_read_mappings" --contig_length_column "contig_length" --otu_column "full_read_mappings_strain" --read_length 100 -o $(RESULTS_PATH)$@ 
	ln -s $(RESULTS_PATH)$@ $@

read_mappings_with_coverage_sorted.tsv: calculate_genome_coverage.py
	time ./$< $(RESULTS_PATH)"read_mappings.csv" --read_mapping_column "full_read_mappings" --contig_length_column "contig_length" --otu_column "full_read_mappings_strain" --read_length 100 --sorted -o $(RESULTS_PATH)$@ 
	ln -s $(RESULTS_PATH)$@ $@
