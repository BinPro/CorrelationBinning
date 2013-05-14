# This experiment aim is to compare how different kmer-lengths affect
# the genomic profile specificity for different taxonomic levels.

DATA_PATH="/proj/b2010008/ProBin/DATA/chris-mock/seperate_assembly/unbalanced"
#DATA_PATH=$HOME"/repos/DATA"

RESULTS_PATH="/proj/b2010008/ProBin/RESULTS/to_big_for_git/2013-03-08/2013-05-14"

KMER_LENGTH=2

score_contigs_multinomial.py:
	ln -s "../../src/programs/score_contigs_multinomial.py" "score_contigs_multinomial.py"

score_mul_0_0_mock_2.tsv: score_contigs_multinomial.py
	time ./score_contigs_multinomial.py $(DATA_PATH)"/contigs.fa" -o $(RESULTS_PATH)"/score_mul_0_0_mock_"$(KMER_LENGTH)".tsv" -t $(DATA_PATH)"/mock_parsed_gen_2_2_complete.txt" -k $(KMER_LENGTH) -d $(DATA_PATH)"/references.fa" --dir_structure "single_fasta_file" --no_taxonomy_info_in_contigs
	ln -s $(RESULTS_PATH)"/score_mul_0_0_mock_"$(KMER_LENGTH)".tsv" "score_mul_0_0_mock_"$(KMER_LENGTH)".tsv"