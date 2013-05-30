## ##
## Author: Brynjar Smari Bjarnason 
## e-mail: binni@binnisb.com

CONTIG_FILES="/proj/b2010008/nobackup/projects/chris-mock/Project_ID793_dAmore/Sample_50ng_even/metassemble/assemblies/velvet/noscaf/newbler/ma-merge.fa"
OUTPUT="/proj/b2010008/ProBin/RESULTS/clustering/2013-04-26_clustering_uppmax"

CONTIG_FILES_MOCK_EVEN_LOW="/proj/b2010008/ProBin/DATA/chris-mock/seperate_assembly/even/even-contigs.fa"


cluster:
	mkdir -p $(OUTPUT)/mock-even-merge-test
	time ProBin.py bin $(CONTIG_FILES) -k 4 -c 41 -a em -o $(OUTPUT)/mock-even-merge-test

mock-even-old:
	mkdir -p $(OUTPUT)/mock-even-low-coverage
	ProBin.py bin $(CONTIG_FILES_MOCK_EVEN_LOW) -k 4 -a em -c 59 -o $(OUTPUT)/mock-even-low-coverage

