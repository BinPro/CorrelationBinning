## ##
## Author: Brynjar Smari Bjarnason 
## e-mail: binni@binnisb.com

CONTIG_FILES="/proj/b2010008/nobackup/projects/chris-mock/Project_ID793_dAmore/Sample_50ng_even/metassemble/assemblies/velvet/noscaf/newbler/ma-merge.fa"
OUTPUT="/proj/b2010008/ProBin/RESULTS/clustering/2013-04-26_clustering_uppmax/mock-even-merge-test"

cluster:
	mkdir -p $(OUTPUT)
	time ProBin.py bin $(CONTIG_FILES) -k 4 -c 41 -a em -o $(OUTPUT)

