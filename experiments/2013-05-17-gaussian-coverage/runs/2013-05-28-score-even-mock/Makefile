DATA_PATH=$(HOME)"/repos/Correlation-Binning/src/tests/fixtures/"

contigs_even_test.fna:
	ln -s $(DATA_PATH)$@ $@

contigs_with_taxonomic_info_test.tsv:
	ln -s $(DATA_PATH)$@ $@

time_series_2013-05-26_12.tsv:
	ln -s $(DATA_PATH)$@ $@

parsed_gen_0_0_mock_test_complete_2013-05-24.txt:
	ln -s $(DATA_PATH)$@ $@

mock_references_test4.fa:
	ln -s $(DATA_PATH)$@ $@

score_time_series.py:
	ln -s "../../src/programs/"$@ $@

score_time_series_test.tsv: score_time_series.py contigs_even_test.fna contigs_with_taxonomic_info_test.tsv time_series_2013-05-26_12.tsv parsed_gen_0_0_mock_test_complete_2013-05-24.txt mock_references_test4.fa
	time ./$< contigs_even_test.fna contigs_with_taxonomic_info_test.tsv time_series_2013-05-26_12.tsv -o $@ -t parsed_gen_0_0_mock_test_complete_2013-05-24.txt -k 2 -d mock_references_test4.fa --first_data "2012-03-25" --last_data "2013-01-18"

