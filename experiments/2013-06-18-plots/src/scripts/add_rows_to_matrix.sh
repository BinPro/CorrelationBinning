#!/usr/bin/env bash
## PARAMETERS:
## -----------
#  $1 - head directory path
#  $2 - level
#  $3 - basename for precision matrix
#  $4 - kmer_length min
#  $5 - kmer_length max
#
## -----------

FILES=$1/*$2
#echo "Add newline to each row..."
#for f in $FILES
#do
#    echo >> $f
#done

echo "Concatenate all lines into a single matrix..."
KMER_MIN=$4
KMER_MAX=$5
KMER_LENGTH=$KMER_MIN
rm $3_$2.tsv
echo "kmer_length & 100 & 1000 & 10000" >> $3_$2.tsv
while [[ $KMER_LENGTH -le $KMER_MAX ]]
do
    a=`cat $1/*_$KMER_LENGTH"_"$2`
    echo $KMER_LENGTH" & "$a >> $3_$2.tsv
    ((KMER_LENGTH = KMER_LENGTH +1))
done

