#!/usr/bin/env python

import sys
import optparse
import random
from Bio import Phylo
from io import StringIO 

parser=optparse.OptionParser()
parser.add_option('-t', '--treefile', help='Treefile in nexus from treetime', type='str')
parser.add_option('-c', '--clusters', default= '',help='file with clusters *.clusters or file with list of nodes', type='str')
parser.add_option('-o', '--outfile', default= '',help='output file with mutations in clusters', type='str')
parser.add_option('-s', '--singletons', default= '',help='file with set of nodes for which I want to get mutations', type='str')
##get options
options, args=parser.parse_args()

tree_with_muts = ""
with open(options.treefile, 'r') as intree:
	for l in intree:
		l = l.rstrip("\n")
		if l.startswith(" Tree tree1="):
			tree_with_muts = l.replace(" Tree tree1=","")
			
#Work with tree
handle = StringIO(tree_with_muts)
tree = Phylo.read(handle, "newick")
#print(tree)

			
			
#read clusters
clusters={} 
if options.clusters != '':
	with open(options.clusters, 'r') as inclust:
		for l in inclust:
			l=l.rstrip("\n")
			if not l=="":
				parts=l.split("\t")
				if not parts[0] in clusters:
					clusters[parts[0]]= list()
				clusters[parts[0]].append(parts[1])

#print(clusters)

#read singletons
singletons  = list()
if options.clusters == '' and options.singletons != '': #I parse clusters if they are provided and parse singletons when there are no clusters
	with open(options.singletons, 'r') as insin:
		for l in insin:
			l=l.rstrip("\n")
			if not l=="":
				parts=l.split("\t")
				if not parts[0] in singletons:
					singletons.append(parts[0])



###work with clusters
with open(options.outfile, 'w') as ouf, open(options.outfile+".vcf", 'w') as oufR:
	oufR.write("##fileformat=VCFv4.1\n##FILTER=<ID=PASS,Description=\"All filters passed\">###\n#CHROM	POS	ID	REF	ALT	QUAL	FILTER 	INFO\n")
	if options.clusters != '':
		for cl in clusters:
			common_ancestor = tree.common_ancestor(clusters[cl])
			Phylo.draw_ascii(common_ancestor)
			#nodes = common_ancestor.get_nonterminals()
			for node in common_ancestor.find_clades():
				if node.comment:
					mutline=node.comment.replace('&mutations="','').replace('"','')
					muts = mutline.split(",")
					ouf.write("%s %s %s\n" %(cl, node.name, mutline))
					for m in muts:
						if not "-" in m:
							if not m=="":
								change=m[:1]+">"+m[-1]
								position=m[1:-1]
								oufR.write("NC_045512.2	%s	.	%s	%s	68	PASS	Cluster=%s,Node=%s\n" %(position, m[:1], m[-1], cl, node.name))
				else:
					ouf.write("%s %s\n" %(cl, node.name))
	#figure out what to do with singletons
	elif options.clusters == '' and options.singletons !='':
		for clade in tree.find_clades():
			if clade.name in singletons:
				if clade.comment:
					mutline=clade.comment.replace('&mutations="','').replace('"','')
					muts = mutline.split(",")
					ouf.write("singleton %s %s\n" %(clade.name, mutline))
					for m in muts:
						if not "-" in m:
							if not m=="":
								change=m[:1]+">"+m[-1]
								position=m[1:-1]
								oufR.write("NC_045512.2	%s	.	%s	%s	68	PASS	Cluster=singleton,Node=%s\n" %(position, m[:1], m[-1],  clade.name))
				else:
					ouf.write("singleton %s\n" %(clade.name))
	elif options.clusters == '' and options.singletons =='':
		for node in tree.find_clades():
			if node.comment:
				mutline=node.comment.replace('&mutations="','').replace('"','')
				muts = mutline.split(",")
				ouf.write("%s %s\n" %(node.name, mutline))
				for m in muts:
					if not "-" in m:
						if not m=="":
							change=m[:1]+">"+m[-1]
							position=m[1:-1]
							oufR.write("NC_045512.2	%s	.	%s	%s	68	PASS	Cluster=All,Node=%s\n" %(position, m[:1], m[-1],  node.name))

