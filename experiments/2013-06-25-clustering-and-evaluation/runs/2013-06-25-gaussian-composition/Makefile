RESULTS_PATH="~/repos/RESULTS/"

even_contigs_good_readmapping_low_coverage.fna.kmer4.signatures:
	ln -s $(HOME)"/Dropbox/ProBinDATA/mock/composition/signatures/"$@

even_contigs_good_readmapping_low_coverage.fna:
	ln -s $(HOME)"/Dropbox/ProBinDATA/mock/composition/fasta/"$@

cluster_kmeans: even_contigs_good_readmapping_low_coverage.fna
	ProBin.py bin composition -a kmeans -c 40 -r 6 -i 20 -k 4 -o $(RESULTS_PATH)"clustering/2013-06-24/" -m isotropic_gaussian $<
	touch $@

cluster_em:  even_contigs_good_readmapping_low_coverage.fna
	ProBin.py bin composition -a em -c 40 -r 6 -i 20 -k 4 -o $(RESULTS_PATH)"clustering/2013-06-24/" -m isotropic_gaussian $<
	touch $@