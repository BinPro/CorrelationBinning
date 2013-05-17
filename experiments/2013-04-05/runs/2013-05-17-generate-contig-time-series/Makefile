# This experiment aim is to trim the kalmar time-series file
# to the number of genomes present in the mock community

RESULTS_PATH=$(HOME)"/repos/RESULTS/to_big_for_git/2013-05-02/low_coverage_sample/even/"

TIME_SERIES_FILE=$(HOME)"/repos/RESULTS/to_big_for_git/time_series/time_series_2013-05-16_59.tsv"

match_contigs_and_time_series.py:
	ln -s "../../src/programs/"$@ $@

contigs_with_samples.tsv: match_contigs_and_time_series.py
	time ./$< $(TIME_SERIES_FILE) --contigs_file $(RESULTS_PATH)"read_mappings_with_coverage_sorted_no_NA.tsv" --first_data "2012-03-25" --last_data "2013-01-18" -o $(RESULTS_PATH)$@ 
	ln -s $(RESULTS_PATH)$@ $@
