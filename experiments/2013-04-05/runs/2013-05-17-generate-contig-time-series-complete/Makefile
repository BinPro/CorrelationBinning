# This experiment aim is to count the reads per 
# unambiguous genome, sort and remove N/A
# values and simulate time series data for each contig 
# based on the kalmar time-series and the within-genome
# distribution of the reads from the even mock.

RESULTS_PATH=$(HOME)"/repos/RESULTS/to_big_for_git/2013-05-02/low_coverage_sample/even/"

TIME_SERIES_FILE=$(HOME)"/repos/RESULTS/to_big_for_git/time_series/time_series_2013-05-16_59.tsv"

calculate_genome_coverage.py:
	ln -s "../../src/programs/"$@ $@

match_contigs_and_time_series.py:
	ln -s "../../src/programs/"$@ $@

read_mappings_with_coverage_sorted_no_NA.tsv: calculate_genome_coverage.py
	time ./$< $(RESULTS_PATH)"read_mappings.csv" --read_mapping_column "full_read_mappings" --contig_length_column "contig_length" --otu_column "full_read_mappings_strain" --read_length 100 --sorted --remove_na -o $(RESULTS_PATH)$@ 
	ln -s $(RESULTS_PATH)$@ $@

contigs_with_samples.tsv: match_contigs_and_time_series.py read_mappings_with_coverage_sorted_no_NA.tsv
	time ./$< $(TIME_SERIES_FILE) --contigs_file read_mappings_with_coverage_sorted_no_NA.tsv --first_data "2012-03-25" --last_data "2013-01-18" -o $(RESULTS_PATH)$@ 
	ln -s $(RESULTS_PATH)$@ $@
