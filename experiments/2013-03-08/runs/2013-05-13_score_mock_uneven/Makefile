# This experiment aim is to compare how different kmer-lengths affect
# the genomic profile specificity for different taxonomic levels.

DATA_PATH="/proj/b2010008/ProBin/DATA/chris-mock/seperate_assembly/unbalanced/"
#DATA_PATH=$HOME"/repos/DATA"

RESULTS_PATH="/proj/b2010008/ProBin/RESULTS/to_big_for_git/2013-03-08/2013-05-14/"

score_contigs_multinomial.py:
	ln -s "../../src/programs/"$@ $@

score_mul_0_0_mock_%.tsv: score_contigs_multinomial.py
	time ./score_contigs_multinomial.py $(DATA_PATH)"contigs.fa" -o $(RESULTS_PATH)$@ -t $(DATA_PATH)"parsed_gen_0_0_mock_complete_fixed.txt" -k $* -d $(DATA_PATH)"references.fa" --dir_structure "single_fasta_file" --no_taxonomy_info_in_contigs
	ln -s $(RESULTS_PATH)$@ $@
