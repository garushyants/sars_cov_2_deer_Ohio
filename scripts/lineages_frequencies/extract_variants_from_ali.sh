#!/bin/bash

#Delta
../../scripts/lineages_frequencies/get_variants_from_ali.py -a Delta.human.deer.reference.mod.fasta -o Delta.variants_from_ali -s
cut -f2 Delta.clusters.final.csv | cut -f1 -d'|' | sed 's/hCoV-19\///' > tmp
grep -f tmp Delta.variants_from_ali.strict.tsv > Delta.variants_from_ali.strict.only_clusters.tsv
grep "deer/USA/OH-OSU" Delta.variants_from_ali.strict.tsv | grep -v -f tmp > Delta.variants_from_ali.strict.only_singletons.tsv
rm tmp
#Alpha
../../scripts/lineages_frequencies/get_variants_from_ali.py -a Alpha.Deer.reference.corrected.fasta -o Alpha.variants_from_ali -s
cut -f2 Alpha.clusters | cut -f1 -d'|' | sed 's/hCoV-19\///' > tmp
grep -f tmp Alpha.variants_from_ali.strict.tsv > Alpha.variants_from_ali.strict.only_clusters.tsv
rm tmp

