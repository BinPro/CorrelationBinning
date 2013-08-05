#!/usr/bin/env bash
## PARAMETERS:
## -----------
#  $1 - head directory path
#  $2 - level
#  $3 - basename for precision matrix
#
## -----------

FILES=$1/*$2
#echo "Add newline to each row..."
#for f in $FILES
#do
#    echo >> $f
#done

echo "Concatenate all lines into a single matrix..."
rm $3_$2.tsv
echo "Number_of_samples & Precision" >> $3_$2.tsv
for N_SAMPLES in 2 5 10 15 30 46
do
    a=`cat $1/*_$N_SAMPLES"_2_"$2`
    echo $N_SAMPLES" & "$a >> $3_$2.tsv
done

