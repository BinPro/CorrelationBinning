# Calculate number of reads mapping to each contig
# given the bedtools coverage.2
DATA_PATH="/proj/b2010008/ProBin/DATA/chris-mock/seperate_assembly/even/original/map/"
RESULTS_PATH="/proj/b2010008/ProBin/DATA/chris-mock/seperate_assembly/even/"
pair_ref-s-md-s.coverage.2:
	ln -s $(DATA_PATH)$@ $@

calculate_nr_reads.py:
	ln -s "../../src/programs/"$@ $@

calculate_genome_coverage.py:
	ln -s "../../src/programs/"$@ $@

contig-purity.tsv:
	ln -s $(RESULTS_PATH)"masmvali/masmvaliout/"$@ $@

reads_per_contig.tsv: calculate_nr_reads.py pair_ref-s-md-s.coverage.2
	./$< pair_ref-s-md-s.coverage.2 -o $@ --read_length 100

match_purity_rpc.py:
	ln -s "../../src/programs/"$@ $@

contig-purity-with-total.tsv:contig-purity.tsv reads_per_contig.tsv match_purity_rpc.py
	./match_purity_rpc.py contig-purity.tsv --rpc reads_per_contig.tsv -o $@