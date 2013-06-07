# This experiment aim is to compare how different kmer-lengths affect
# the genomic profile specificity for different taxonomic levels.

DATA_PATH="/proj/b2010008/ProBin/DATA/chris-mock/seperate_assembly/even/"
DATA_PATH2="/proj/b2010008/ProBin/RESULTS/to_big_for_git/2013-05-02/low_coverage_sample/even/"
#DATA_PATH=$HOME"/repos/DATA"

RESULTS_PATH="/proj/b2010008/ProBin/RESULTS/to_big_for_git/2013-03-08/2013-06-03/"

score_contigs_multinomial.py:
	ln -s "../../src/programs/"$@ $@

contigs.fa:
	ln -s $(DATA_PATH)"even-contigs.fa" $@

parsed_gen_0_0_mock_complete_fixed3.txt:
	ln -s $(DATA_PATH2)$@ $@

references.fa:
	ln -s $(DATA_PATH)$@ $@

contigs_with_taxonomic_info.tsv:
	ln -s $(DATA_PATH2)$@ $@

score_mul_0_0_mock_%.tsv: score_contigs_multinomial.py contigs.fa parsed_gen_0_0_mock_complete_fixed3.txt references.fa contigs_with_taxonomic_info.tsv
	time ./score_contigs_multinomial.py contigs.fa -o $(RESULTS_PATH)$@ -t parsed_gen_0_0_mock_complete_fixed3.txt -k $* -d references.fa --dir_structure "single_fasta_file" --no_taxonomy_info_in_contigs --filter_file contigs_with_taxonomic_info.tsv
	ln -s $(RESULTS_PATH)$@ $@
