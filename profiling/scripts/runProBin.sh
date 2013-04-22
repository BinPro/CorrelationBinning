#/usr/bin/env bash
DIREC=~/MasterProject/ProBin
kernprof.py --line-by-line $DIREC/src/ProBin.py $DIREC/tests/generated_contigs_10000_test.fna -c 4 -mc dirichlet
#kernprof.py --line-by-line $DIREC/src/ProBin.py $DIREC/tests/generated_contigs_test.fna -c 4 -mc dirichlet
