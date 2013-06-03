#!/bin/bash
for f in $@
do
  echo Cluster `echo $f | sed -re "s/.*(00..).fa/\1/"`:$(sed "s/>/,/g" $f ) 
done

