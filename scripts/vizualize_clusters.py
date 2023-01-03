#!/usr/bin/env python

import sys
import ete3
import optparse
from ete3 import Tree, NodeStyle, AttrFace, faces, TreeStyle, SeqMotifFace,TextFace
import random
import re

parser=optparse.OptionParser()
parser.add_option('-i', '--infile', help='Treefile to visualize', type='str')
parser.add_option('-m', '--mutations', default = '', help ='vcf file with annotated mutations',type='str')
#parser.add_option('-t', '--nonbinary', action="store_true", default=False)
parser.add_option('-c', '--clusterIDs', help='file with cluster descriptions',default = '', type = 'str')
parser.add_option('-r', '--root', default= "", help='', type='str')
parser.add_option('-f','--format',default = 1, help='tree fromat for ete3, default =0',type = 'int')
parser.add_option('-o', '--outfolder', help='path to folder where to store pdfs', default = 'cluster_figures', type='str')

##get options
options, args=parser.parse_args()
tree = ete3.Tree(options.infile, format=options.format) #input tree is from treetime and nonbinary
outfolder = options.outfolder
#nbchecker = options.nonbinary
clusterIDs=options.clusterIDs
mutationsfile=options.mutations
root=options.root

aminoacidscodefiles= "../../data/genetic_code.txt"

#cutoff is needed if tree is binary
cutoff = 1e-05

if not root == "":
	tree.set_outgroup(tree&root)

############################
#read file with aminoacids
aminoacidscode =dict()
with open(aminoacidscodefiles, 'r') as aafile:
	for l in aafile:
		l=l.rstrip('\n')
		codon,one,three=l.split()
		aminoacidscode[three]=one

print(aminoacidscode)

#read file with mutations
MutDict = dict()
with open(mutationsfile, 'r') as mutf:
	#next(mutf)
	for ll in mutf:
		if not ll.startswith('#'):
			ll=ll.rstrip('\n')
			main = ll.split('\t')[7]
			print(main)
			parts = re.split(',|;|\|', main)
			print(parts)		
			aachange = re.split('/|\.',parts[5])[1]
			gene = parts[7]
			cluster =parts[0].replace('Cluster=','')
			node = parts[1].replace('Node=','').replace('\\','|')
			print(aachange, gene, cluster, node, sep =' ')
			if parts[2].startswith('EFF=missense_variant') or parts[2].startswith('EFF=stop_gained'):
				aachangefinal = aminoacidscode[aachange[0:3]]+aachange[3:-3]
				if not aachange.endswith('*'):
					aachangefinal += aminoacidscode[aachange[-3:]]
				else:
					aachangefinal += aachange[5:]
				if not cluster in MutDict:
					MutDict[cluster] =dict()
				if not node in MutDict[cluster]:
					MutDict[cluster][node] = dict()
				if not gene in MutDict[cluster][node]:
					MutDict[cluster][node][gene] = list()
				MutDict[cluster][node][gene].append(aachangefinal)
#print(MutDict)	
############################
####Get clusters from file and get the correct clusternodes
clusters=dict()
with open(clusterIDs, 'r') as clusterfile:
	for l in clusterfile:
		l=l.rstrip("\n")
		if not l == '':
			cluster,leave=l.split("\t")
			if not cluster in clusters:
				clusters[cluster]=list()
			clusters[cluster].append(tree.get_leaves_by_name(leave)[0])
	
##Visualize
def layout(node):
	if node.is_leaf():
		if "/deer/" in node.name:
			deerdesc = faces.AttrFace("name", fsize=8)
			deerdesc.margin_left = 25
			faces.add_face_to_node(deerdesc, node, 0, aligned=True)
		else:
			if "_" in node.name:
				node.name = node.name.replace("_"," ")
			humandesc = faces.AttrFace("name", fsize=7)
			humandesc.margin_left = 5
			humandesc.margin_bottom = 1
			humandesc.margin_top = 1
			humandesc.margin_right = 5
			faces.add_face_to_node(humandesc, node, 0, aligned = False)

ts = TreeStyle()
#ts.mode = "c"
ts.scale = 1000000
#ts.tree_width = 500
ts.show_leaf_name = False
#ts.draw_guiding_lines =True
#ts.guiding_lines_type = 2
#ts.show_branch_support = True
ts.layout_fn = layout

widthsingl = 2
widthcluster = 3

colorlist = ["#8dd3c7","#ffffb3","#bebada","#fb8072","#80b1d3","#fdb462","#b3de69","#fccde5","#d9d9d9","#bc80bd","#ccebc5","#ffed6f","#7fc97f","#beaed4","#fdc086","#ffff99","#386cb0","#f0027f"]

j=0
for clid in clusters:
	commonnode = tree.get_common_ancestor(clusters[clid])
	print(commonnode)
	MutClDict = MutDict[clid]
	print(MutClDict)
	for node in MutClDict:
		#text = ""
		branch = commonnode.search_nodes(name=node)[0]
		for gene in MutClDict[node]:
			text = gene+":"+ ",".join(MutClDict[node][gene])
			exec('branch.add_face(TextFace("'+text+'",fsize=8), column =0, position="branch-top")')
		#text = text.rstrip(";")
		#print(branch)
		#exec('branch.add_face(TextFace("'+text+'",fsize=8),column =1, position="branch-top")')
	#set styles for clusters
	stylename = "nsStyle"+str(clid)
	leafstylename ="lsStyle"+str(clid)
	exec(stylename+"=NodeStyle()")
	color = colorlist[j]+ "66" #"%06x" % random.randint(0, 0xFFFFFF) + "66" #Additional digits are for transparency
	#print("%i %s" %(clid,color))
	#exec(stylename+"['bgcolor'] = '"+color+"'")
	exec(stylename+"['hz_line_color'] = '"+color+"'")
	exec(stylename+"['vt_line_color'] = '"+color+"'")
	exec(stylename+"['vt_line_width'] = "+str(widthcluster))
	exec(stylename+"['hz_line_width'] = "+str(widthcluster))
	exec(stylename+"['size'] = 0")
	
	exec(leafstylename+"=NodeStyle()")
	exec(leafstylename+"['vt_line_width'] = "+str(widthcluster))
	exec(leafstylename+"['hz_line_width'] = "+str(widthcluster))
	exec(leafstylename+"['hz_line_color'] = '"+color+"'")
	exec(leafstylename+"['vt_line_color'] = '"+color+"'")
	exec(leafstylename+"['size'] = 0")
	
	
	###use it to set styles for clusters
	for n in commonnode.traverse():
		exec("n.set_style(lsStyle"+str(clid)+")")
	exec("commonnode.set_style(nsStyle"+str(clid)+")")
	#commonnode.show()
	j+=1
	commonnode.render(outfolder+"/"+clid+".svg", w = 1000, dpi =300, units= 'px', tree_style = ts)




