DATA_PATH="/proj/b2010008/ProBin/DATA/chris-mock/"
DATA_PATH2="/proj/b2010008/ProBin/RESULTS/to_big_for_git/2013-05-02/low_coverage_sample/even/cross_validation_and_real/cv/"
DATA_PATH3="/proj/b2010008/ProBin/RESULTS/to_big_for_git/"
RESULTS_PATH="/proj/b2010008/ProBin/RESULTS/to_big_for_git/2013-06-04/2013-06-04/"

even-contigs.fa:
	ln -s $(DATA_PATH)"seperate_assembly/even/"$@ $@

contigs_with_taxonomic_info.tsv:
	ln -s $(DATA_PATH2)$@ $@

time_series_2013-06-04_46-cv.tsv:
	ln -s $(DATA_PATH3)"time_series/"$@ $@

parsed_gen_0_0_mock_complete_fixed3_6.txt:
	ln -s $(DATA_PATH2)$@ $@

references.fa:
	ln -s $(DATA_PATH)"references/"$@ $@

score_combined.py:
	ln -s "../../src/programs/"$@ $@

score_combined_cv_46_columns_6_10000.tsv: score_combined.py even-contigs.fa contigs_with_taxonomic_info.tsv time_series_2013-06-04_46-cv.tsv parsed_gen_0_0_mock_complete_fixed3_6.txt references.fa
	time ./$< even-contigs.fa contigs_with_taxonomic_info.tsv time_series_2013-06-04_46-cv.tsv -o $(RESULTS_PATH)$@ -t parsed_gen_0_0_mock_complete_fixed3_6.txt -k 6 -f 10000 -d references.fa --first_data "2012-03-25" --last_data "2013-01-18"
	ln -s $(RESULTS_PATH)$@ $@