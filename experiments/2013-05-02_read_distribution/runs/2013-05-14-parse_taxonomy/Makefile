## ##
## Author: Johannes Alneberg 
## e-mail: johannesalneberg@gmail.com
## This makefile executes the generation of 
## taxonomy data for the mock community data.

PHYLOGENY_FILE="/proj/b2010008/ProBin/DATA/chris-mock/references/phylogeny-references.tsv"

parse_taxonomy_complete.py:
	ln -s ../../../2013-01-18/src/programs/parse_taxonomy_complete.py

parsed_gen_0_0_mock_complete.txt: parse_taxonomy_complete.py
	./parse_taxonomy_complete.py $(PHYLOGENY_FILE) -o parsed_gen_0_0_mock -s 0 -g 0 --file_name_column "fasta_name" --species_column "topname" --backup_file_name_column "xls_name_in_info_run" --backup_species_column "xls_name_in_info_run"