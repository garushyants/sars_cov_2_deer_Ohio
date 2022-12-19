#!/bin/bash

NUMBERMUT=$1
CYCLES=$2
FOLDER=subsamples_human_${NUMBERMUT}_${CYCLES}

grep -v "#" UShER_all.vcf > tmp_file
grep "#" UShER_all.vcf > tmp_header

mkdir $FOLDER
cd $FOLDER

for (( c=1; c<=$CYCLES; c++ ))
do
echo $c
shuf -n $NUMBERMUT ../tmp_file > all${c}.vcf 
cat ../tmp_header all${c}.vcf > all${c}.header.vcf 
rm all${c}.vcf
done

cd ..
rm tmp_file
rm tmp_header


