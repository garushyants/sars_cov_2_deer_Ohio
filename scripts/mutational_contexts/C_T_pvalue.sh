#!/bin/bash

NUMBERMUT=435
CYCLES=1000
REALDATA='./Delta.clusters.snpeff.nopr.forMutationalSignatures.vcf'

FILE=calculate_pvalue_${NUMBERMUT}_${CYCLES}.txt

grep -v "#" UShER_all.vcf > tmp_file

cat > ${FILE}
for (( c=1; c<=$CYCLES; c++ ))
do
echo $c
shuf -n $NUMBERMUT tmp_file | grep -P "\tC\tT\t" | wc -l >> $FILE

done

rm tmp_file
###
REALVALUE=$(grep -v "#" $REALDATA | grep -P "\tC\tT\t" | wc -l)
echo $REALVALUE
###
awk '{ sum += $1; n++ } END { if (n > 0) print sum / n; }' $FILE

Rscript -e 'args <- commandArgs(TRUE)' \
   -e 'perc.rank <- function(x, xo)  1-length(x[x < xo])/length(x)' \
	-e 'd<-scan(args[2], quiet=TRUE)' \
      -e 'perc.rank(d,args[1])' ${REALVALUE} ${FILE}


