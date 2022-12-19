#!/usr/bin/env python

import sys
import optparse
import random
import re

#The goal of this script is to correct snpeff problems for situations when ref does not match genome
parser=optparse.OptionParser()
parser.add_option('-v', '--vcf', default = '', help ='input vcf file annotated with SnpEff. The recomended parameters are: -v -canon -formatEff -no-upstream -no-downstream',type='str')
#parser.add_option('-o', '--outvcf', help='Name of ouput vcf file. No default is provided', default = '', type='str')
parser.add_option('-c', '--geneticCode', help='Path to file with genetic code', default='/panfs/pan1.be-md.ncbi.nlm.nih.gov/virbac/sars_cov_2_deers/scripts/genetic_code.txt', type='str')

##get options
options, args=parser.parse_args()
vcf = options.vcf
geneticCodeFile=options.geneticCode
############################
#read genetic code
GeneCodeDict =dict()
with open(geneticCodeFile, 'r') as gencin:
	for l in gencin:
		l=l.rstrip('\n')
		codon,sh,lo = l.split(' ')
		GeneCodeDict[codon]=lo


############################

outvcflist = list()
"method to read and parse vcf files"
with open(vcf, 'r') as mutf:
	for l in mutf:
		l=l.rstrip('\n')
		if "WARNING_REF_DOES_NOT_MATCH_GENOME" in l:
			print(l)
			parts = l.split('\t')
			position=parts[1]
			refbase=parts[3]
			altbase=parts[4]
			annotation=parts[7]
			preparts=annotation.replace(';',',').split(',')
			annot_parts = preparts[2].split('|')
			muttype=annot_parts[0].replace("EFF=","").split('(')[0]
			gene = annot_parts[5]
			prcoding = annot_parts[6]
			#
			change = annot_parts[2]
			refcodon,altcodon=change.split('/')
			#
			###This part is here to correct warnings for warnings
			if len(annot_parts)>11 and prcoding == 'protein_coding':
				if 'WARNING_REF_DOES_NOT_MATCH_GENOME' in annot_parts[11]:
					correctrefcodon=re.sub(r'[ATGC]', refbase, refcodon)
					if GeneCodeDict[correctrefcodon.upper()] == GeneCodeDict[altcodon.upper()]:
						muttype='EFF=synonymous_variant(LOW|SILENT'
					else:
						if GeneCodeDict[correctrefcodon.upper()] == 'Ter':
							muttype = 'EFF=stop_lost(HIGH|NONSENSE'
						elif GeneCodeDict[altcodon.upper()] == 'Ter':
							muttype = 'EFF=stop_gained(HIGH|NONSENSE'
						else:
							muttype='EFF=missense_variant(MODERATE|MISSENSE'
					protchange=annot_parts[3].replace(GeneCodeDict[refcodon.upper()],GeneCodeDict[correctrefcodon.upper()],1)
			correctedline="\t".join(parts[0:7])+"\t"+preparts[0]+','+preparts[1]+";"+muttype+"|"+correctrefcodon+"/"+altcodon+"|"+protchange+"|"+"|".join(annot_parts[4:11])
			outvcflist.append(correctedline)
							
		else:
			outvcflist.append(l)
			
			
#write to file
outfile=vcf.replace(".vcf",".corrected.vcf")
#print(outfile)
with open(outfile, 'w') as ouf:
	for li in outvcflist:
		ouf.write("%s\n" % li)
	print("Corrected vcf printed to "+outfile)
			
#print(ouvcf)
	
	



