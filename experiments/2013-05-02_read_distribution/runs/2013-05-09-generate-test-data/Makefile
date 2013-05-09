## ##
## Author: Johannes Alneberg 
## e-mail: johannesalneberg@gmail.com
## This makefile executes the generation of 
## test data for the mock community data.

PHYLOGENY_FILE="/proj/b2010008/ProBin/DATA/chris-mock/test_data/phylogeny_mock_test.tsv"

parsed_gen_0_0_mock_test_complete.txt:
	./parse_taxonomy_complete.py $(PHYLOGENY_FILE) -o parsed_gen_0_0_mock_test -s 0 -g 0 --file_name_column "fasta_name" --species_column "topname"