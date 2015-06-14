#!/bin/bash

# Systems Immunology Lab, WUSTL

NAME=$1  ## common prefix in BED files to be processed, e.g. [proB]_H3K27ac, [proB]_H3K4me1, etc
WFILE=$2 ## genomic windows
CHROMSIZES=$3

echo "This script is intended to run on multiple CPUs"

KK=`awk '{print $1}' $CHROMSIZES`

for i in $NAME*.bed
do
  bedtools intersect -c -f 0.5 -b $i -a hg19.200bp_windows.bed > ${i%%.bed}_binary.bed &
done
wait

for i in $NAME*_binary.bed
do
  for j in $KK
  do
    grep -P "$j\t" $i | awk '{print $4}' > ${i%%_binary.bed}_${j}.x &
  done
  wait
  echo "done processing $i"
done

ls $NAME*_chr1.x | sed "s/_chr1\.x//g" |  sed "s/${NAME}_//g"> $NAME.$$.marks

PP=`cat $NAME.$$.marks`

for j in $KK
do
  echo "$NAME $j" | sed "s/ /\t/g" > ${NAME}_${j}_binary.txt
  echo $PP | sed "s/ /\t/g" >>  ${NAME}_${j}_binary.txt
  paste $NAME*_${j}.x >> ${NAME}_${j}_binary.txt
done

rm *.x *binary.bed $NAME.$$.marks

