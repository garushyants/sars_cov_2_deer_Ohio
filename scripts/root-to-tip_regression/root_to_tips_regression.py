#!/usr/bin/env python

import os
import sys
from operator import add
import ete3
import optparse
import re


parser=optparse.OptionParser()
parser.add_option('-t', '--tree', help='Treefile in newick format', type='str')
parser.add_option('-v', '--vcf', help='Input annotated vcf file', type='str')
parser.add_option('-f','--informat',default = 1, help='tree format for input tree: default format is 0 format in ete3',type = 'int')
parser.add_option('-o', '--outfile', help='output file with mutation counts', default = '', type='str')
parser.add_option('-g', '--gene', help='Gene to count mutations for', default = 'S', type='str')
parser.add_option('-w', '--correctWarnings', help='This option allows to correct for WARNING_REF_DOES_NOT_MATCH_GENOME', action="store_true", default=False)
parser.add_option('-c', '--geneticCode', help='Path to file with genetic code', default='../data/genetic_code.txt', type='str')

##get options
options, args=parser.parse_args()
tree = ete3.Tree(options.tree, format=options.informat) #input tree is from treetime and nonbinary
outfile = options.outfile
vcffile=options.vcf
genesymb = options.gene
checkWarnings=options.correctWarnings
geneticCodeFile=options.geneticCode

GeneCodeDict =dict()
if checkWarnings:
	with open(geneticCodeFile, 'r') as gencin:
		for l in gencin:
			l=l.rstrip('\n')
			codon,sh,lo = l.split(' ')
			GeneCodeDict[codon]=lo
		


VCFDict = dict()
#Read vcf
with open(vcffile, 'r') as vcfin:
	for l in vcfin:
		if not l.startswith("#"):
			l=l.rstrip('\n')
			annotation = l.split('\t')[7]
			preparts=annotation.replace(';',',').split(',')
			node=preparts[1].replace('Node=','')
			annot_parts = preparts[2].split('|')
			muttype=annot_parts[0].replace("EFF=","").split('(')[0]
			gene = annot_parts[5]
			prcoding = annot_parts[6]
			###This is an added part to correct for warnings
			if checkWarnings:
				if len(annot_parts)>11 and prcoding == 'protein_coding':
					if 'WARNING_REF_DOES_NOT_MATCH_GENOME' in annot_parts[11]:
						change = annot_parts[2]
						refcodon,altcodon=change.split('/')
						altcodon=altcodon.upper()
						refbase=l.split('\t')[3]
						correctrefcodon=re.sub(r'[ATGC]', refbase, refcodon).upper()
						if GeneCodeDict[correctrefcodon] == GeneCodeDict[altcodon]:
							muttype='synonymous_variant'
						else:
							if GeneCodeDict[correctrefcodon] == 'Ter':
								muttype = 'stop_lost'
							elif GeneCodeDict[altcodon] == 'Ter':
								muttype = 'stop_gained'
							else:
								muttype='missense_variant'
						print('%s\t%s\t%s\t%s' %(refcodon,correctrefcodon,altcodon,muttype))
			####
			if not node in VCFDict:
				VCFDict[node] = dict()
			if not muttype in VCFDict[node]:
				VCFDict[node][muttype]=list()
			VCFDict[node][muttype].append(gene)
			
#print(VCFDict)


###Let's traverse tree and do the counts for each leaf
RegressionDict = dict()

for node in tree.traverse("postorder"):
	if node.name in VCFDict:
		synmutcount = 0
		mismutcount = 0
		Ssyn = 0
		Smis = 0
		allmutcount = 0 
		for muttypes in VCFDict[node.name]:
			if 'synonymous_variant' in VCFDict[node.name]:
				synmutcount = len(VCFDict[node.name]['synonymous_variant'])
				if genesymb in VCFDict[node.name]['synonymous_variant']:
					Ssyn = VCFDict[node.name]['synonymous_variant'].count(genesymb)
			if 'missense_variant' in VCFDict[node.name]:
				mismutcount = len(VCFDict[node.name]['missense_variant'])
				if genesymb in VCFDict[node.name]['missense_variant']:
					Smis = VCFDict[node.name]['missense_variant'].count(genesymb)
			allmutcount += len(VCFDict[node.name][muttypes])
		mutcounts = [allmutcount,synmutcount,mismutcount,Ssyn,Smis]
		for leaf in node:
			if not leaf.name in RegressionDict:
				RegressionDict[leaf.name] = list()
				RegressionDict[leaf.name] = mutcounts
			else:
				mutrecount = list( map(add, RegressionDict[leaf.name], mutcounts) )
				RegressionDict[leaf.name] = mutrecount
				
#print (RegressionDict)

with open(outfile, 'w') as ouf:
	ouf.write('Leaf\tAll\tSyn\tMis\t'+genesymb+'syn\t'+genesymb+'mis\n')
	for k in RegressionDict:
		ouf.write(k +"\t")
		ouf.write("\t".join(map(str, RegressionDict[k])) +"\n")

		
