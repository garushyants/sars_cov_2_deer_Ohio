#!/usr/bin/env python

import sys
import ete3
import optparse
from ete3 import Tree, NodeStyle, AttrFace, faces, TreeStyle, SeqMotifFace,TextFace
import random

parser=optparse.OptionParser()
parser.add_option('-i', '--infile', help='Treefile to visualize', type='str')
parser.add_option('-l', '--leafnames', default= '',help='', type='str')
parser.add_option('-m', '--mutations', default = '', help ='csv file with annotated mutations',type='str')
parser.add_option('-t', '--nonbinary', action="store_true", default=False)
parser.add_option('-c', '--clusterIDs', default = 'deer/USA/OH-OSU', type = 'str')
parser.add_option('-r', '--root', default= "", help='', type='str')
parser.add_option('-f','--format',default = 0, help='default format is 0 format in ete3',type = 'int')
parser.add_option('-o', '--outfolder', help='path to final pdf', default = 'cluster_figures', type='str')
parser.add_option('-n', '--node', help='use internal nodes names as cluster ids', action="store_true", default=False)
parser.add_option('-e', '--onlyReroot', action="store_true", default=False)

##get options
options, args=parser.parse_args()
tree = ete3.Tree(options.infile, format=options.format) #input tree is from treetime and nonbinary
root = options.root
outfolder = options.outfolder
nbchecker = options.nonbinary
clnodeid = options.node
onlyReroot=options.onlyReroot
clusterIDs=options.clusterIDs
mutationsfile=options.mutations

aminoacidscodefiles= "../data/genetic_code.txt"

#cutoff is needed if tree is binary
cutoff = 1e-05

if not root == "":
	tree.set_outgroup(tree&root)
	
if onlyReroot:
	tree.write(format=0, outfile=options.infile+"_rerooted.nwk")
	exit()
	

def convert_to_nonbinary(t, threshold):
	for node in t.iter_descendants("postorder"):
		#print node.dist
		if node.dist < threshold:
			if not node.is_leaf():
				for child in node.children:
					(node.up).add_child(child, dist = child.dist)
				node.detach()
	return(t)

if not nbchecker: #check whether the input tree is nonbinary
	tree = convert_to_nonbinary(tree, cutoff)
	#tree.write(format=1, outfile=options.infile+"_nonbinary.nwk")
	
############################
#read file with aminoacids
aminoacidscode =dict()
with open(aminoacidscodefiles, 'r') as aafile:
	for l in aafile:
		l=l.rstrip('\n')
		codon,three,one=l.split()
		aminoacidscode[three]=one

#read file with mutations
MutDict = dict()
with open(mutationsfile, 'r') as mutf:
	next(mutf)
	for ll in mutf:
		ll=ll.rstrip('\n')
		parts = ll.split(',')
		aachange = parts[7]
		gene = parts[5]
		cluster =parts[0]
		node = parts[2].replace('\\','|')
		aachangefinal = aminoacidscode[aachange[2:5]]
		if not aachange.endswith('*'):
			aachangefinal += aachange[5:-3]+aminoacidscode[aachange[-3:]]
		else:
			aachangefinal += aachange[6:]
		if not cluster in MutDict:
			MutDict[cluster] =dict()
		if not node in MutDict[cluster]:
			MutDict[cluster][node] = dict()
		if not gene in MutDict[cluster][node]:
			MutDict[cluster][node][gene] = list()
		MutDict[cluster][node][gene].append(aachangefinal)
#print(MutDict)	
############################

####Go around the tree and find clusters of Ohio deers

clusters = {}
num = 1

for node in tree.iter_descendants("postorder"):
	if node.is_leaf():
		if clusterIDs in node.name:
			node.add_feature('state',1)
		else:
			node.add_feature('state',0)
	else:
		node.state = 1
		ohdeer = []
		nonohdeer = 0
		for no in node.children:
			if no.state == 1:
				ohdeer.append(no)
			else:
				nonohdeer +=1
		if nonohdeer > 0:
			node.state = 0
			if len(ohdeer) > 0:
				for child in ohdeer:
					if not child.is_leaf() and child.dist > cutoff:   ### the child is a cluster founder if it is internal (has at least two terminal descendants) and russian
						if clnodeid:
							clusters[node.name] = child.get_leaf_names()
						else:
							clusters[num] = child.get_leaf_names()
							num +=1
	
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





