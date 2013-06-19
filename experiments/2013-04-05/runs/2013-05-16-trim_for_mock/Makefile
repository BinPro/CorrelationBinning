# This experiment aim is to trim the kalmar time-series file
# to the number of genomes present in the mock community

RESULTS_PATH=$(HOME)"/repos/RESULTS/to_big_for_git/time_series/"

trim_time_series.py:
	ln -s "../../src/programs/"$@ $@

time_series_2013-05-16_59.tsv: trim_time_series.py
	time ./$< $(RESULTS_PATH)"time_series_2013-03-05_184.tsv" --first_data "2012-03-25" --last_data "2013-01-18" -o $(RESULTS_PATH)$@ --number_of_genomes 59
	ln -s $(RESULTS_PATH)$@ $@

time_series_2013-05-16_15.tsv: trim_time_series.py
	time ./$< $(RESULTS_PATH)"time_series_2013-03-05_184.tsv" --first_data "2012-03-25" --last_data "2013-01-18" -o $(RESULTS_PATH)$@ --number_of_genomes 15
	ln -s $(RESULTS_PATH)$@ $@

time_series_2013-05-16_40.tsv: trim_time_series.py
	time ./$< $(RESULTS_PATH)"time_series_2013-03-05_184.tsv" --first_data "2012-03-25" --last_data "2013-01-18" -o $(RESULTS_PATH)$@ --number_of_genomes 40
	ln -s $(RESULTS_PATH)$@ $@

