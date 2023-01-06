#!/usr/bin/env python

import os
import sys
import optparse
from Bio import AlignIO

parser=optparse.OptionParser()
parser.add_option('-a', '--alignment', help='Input alignment in fasta format', type='str')
parser.add_option('-o', '--outprefix', default= '',help='output files prefix', type='str')
parser.add_option('-r', '--reference', default = 'NC_045512.2', help='ID of reference strain before the first space', type ='str')
parser.add_option('-s', '--short', action="store_true", default=False)

options, args=parser.parse_args()
alignment = AlignIO.read(options.alignment, "fasta")
outprefix=options.outprefix
ref=options.reference

AlignDict=dict()

counter = 0
ref_num = -5
for record in alignment:
	seqname = ""
	if options.short:
		name_parts = str(record.id).split('|')
		seqname=name_parts[0].replace("hCoV-19/","")
		for n in name_parts:
			if "." in n:
				seqname += " "+n
				break
	else:
		seqname = record.id
	AlignDict[counter]=seqname
	if record.id == ref:
		ref_num = counter
	counter +=1
#print(AlignDict)
#print(alignment[ref_num].seq)
with open(outprefix+".all.tsv", 'w') as outall, open(outprefix+".strict.tsv", 'w') as outstrict:
	for j in range(0,len(alignment[ref_num].seq)):
		ref_letter = alignment[ref_num].seq[j]
		#print(ref_letter)
		column = alignment[:,j]
		#print(column)
		pos_in_str = 0

		for l in column:
			if l != '-' and l != 'n' and l != ref_letter:
				#print(pos_in_str)
				#print(l)
				#print(ref_letter)
				outall.write("%i %s %s %s\n" %(j+1, ref_letter.upper(), l.upper(), AlignDict[pos_in_str]))
				if l in ['a','t','g','c']:
					outstrict.write("%i %s %s %s\n" %(j+1, ref_letter.upper(), l.upper(), AlignDict[pos_in_str]))
			pos_in_str +=1
	
