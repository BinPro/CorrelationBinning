DATA_PATH="/proj/b2010008/ProBin/DATA/chris-mock/"
DATA_PATH2="/proj/b2010008/ProBin/RESULTS/to_big_for_git/2013-05-02/low_coverage_sample/even/"
DATA_PATH3="/proj/b2010008/ProBin/RESULTS/to_big_for_git/"
RESULTS_PATH="/proj/b2010008/ProBin/RESULTS/to_big_for_git/2013-05-17/2013-05-28/"

even-contigs.fa:
	ln -s $(DATA_PATH)"seperate_assembly/even/"$@ $@

contigs_with_taxonomic_info.tsv:
	ln -s $(DATA_PATH2)$@ $@

time_series_2013-05-16_59.tsv:
	ln -s $(DATA_PATH3)"time_series/"$@ $@

parsed_gen_0_0_mock_complete_fixed3.txt:
	ln -s $(DATA_PATH)$@ $@

references.fa:
	ln -s $(DATA_PATH)"references/"$@ $@

score_time_series.py:
	ln -s "../../src/programs/"$@ $@

score_time_series_even_only_iso.tsv: score_time_series.py even-contigs.fa contigs_with_taxonomic_info.tsv time_series_2013-05-16_59.tsv parsed_gen_0_0_mock_complete_fixed3.txt references.fa
	time ./$< even-contigs.fa contigs_with_taxonomic_info.tsv time_series_2013-05-16_59.tsv -o $(RESULTS_PATH)$@ -t parsed_gen_0_0_mock_complete_fixed3.txt -k 2 -d references.fa --first_data "2012-03-25" --last_data "2013-01-18"
	ln -S $(RESULTS_PATH)$@ $@