# Executes the reduction of columns in the time-series file.
RESULTS_PATH=$(HOME)"/repos/RESULTS/to_big_for_git/time_series/"

select_columns_rand.py:
	ln -s "../../src/programs/"$@ $@

time_series_2013-05-16_40.tsv:
	ln -s $(RESULTS_PATH)$@ $@

time_series_2013-06-01_15.tsv:
	ln -s $(RESULTS_PATH)$@ $@

all: select_columns_rand.py time_series_2013-05-16_40.tsv
	./$< time_series_2013-05-16_40.tsv -o $(RESULTS_PATH)"time_series_2013-06-01" --first_data "2012-03-25" --last_data "2013-01-18"

a_few_more: select_columns_rand.py time_series_2013-06-01_15.tsv
	./$< time_series_2013-06-01_15.tsv -o $(RESULTS_PATH)"time_series_2013-06-05" --first_data "2012-03-25" --last_data "2013-01-18" --nr_columns 15 10 5 2