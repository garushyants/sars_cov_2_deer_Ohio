#!/usr/bin/env python

import sys
import ete3
import optparse
from ete3 import Tree, NodeStyle, AttrFace, faces, TreeStyle, SeqMotifFace
import random

parser=optparse.OptionParser()
parser.add_option('-i', '--infile', help='Treefile to visualize', type='str')
parser.add_option('-l', '--leafnames', default= '',help='', type='str')
parser.add_option('-t', '--nonbinary', action="store_true", default=False)
parser.add_option('-c', '--clusterIDs', default = 'deer/USA/OH-OSU', type = 'str')
parser.add_option('-r', '--root', default= "", help='', type='str')
parser.add_option('-f','--format',default = 0, help='default format is 0 format in ete3',type = 'int')
parser.add_option('-o', '--outfile', help='path to final pdf', default = 'tree_colored.pdf', type='str')
parser.add_option('-n', '--node', help='use internal nodes names as cluster ids', action="store_true", default=False)
parser.add_option('-e', '--onlyReroot', action="store_true", default=False)
parser.add_option('-u', '--clusteroutfile', default = '', type = 'str')

##get options
options, args=parser.parse_args()
tree = ete3.Tree(options.infile, format=options.format) #input tree is from treetime and nonbinary
root = options.root
outfile = options.outfile
nbchecker = options.nonbinary
clnodeid = options.node
onlyReroot=options.onlyReroot
clusterIDs=options.clusterIDs
clustersoutfile=options.clusteroutfile

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
###	
#Print clusters and save to file
clouf = ''
if clustersoutfile == '':
	clouf = options.infile+".clusters"
else:
	clouf = clustersoutfile
with open(clouf, 'w') as outf:
	for k in clusters:
			####print clusters
		for elem in clusters[k]:
			print (str(k)+"\t"+elem.replace("\'",""))
			outf.write(str(k)+"\t"+elem.replace("\'","")+"\n")
	
##Visualize
def layout(node):
	if node.is_leaf():
		if "/deer/" in node.name:
			deerdesc = faces.AttrFace("name", fsize=10)
			deerdesc.margin_left = 25
			faces.add_face_to_node(deerdesc, node, 0, aligned=False)
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
ts.branch_vertical_margin = 0.1
#ts.root_opening_factor = 1
ts.arc_start = 180 # 0 degrees = 3 o'clock
ts.arc_span = 355
#ts.mode = "c"
#ts.scale = 2000000
ts.tree_width = 3000
ts.show_leaf_name = False
#ts.draw_guiding_lines =True
#ts.guiding_lines_type = 2
#ts.show_branch_support = True
ts.layout_fn = layout

widthsingl = 10 #2
widthcluster = 10 #2

singletoncolor = "#1b9e77" #"%06x" % random.randint(0, 0xFFFFFF)
otherdeercolor = "#d95f02"
ohdeerStyle = NodeStyle()
ohdeerStyle['hz_line_color'] = singletoncolor
ohdeerStyle['vt_line_color'] = singletoncolor
ohdeerStyle['hz_line_width'] = widthsingl
ohdeerStyle['vt_line_width'] = widthsingl
ohdeerStyle["fgcolor"] = "black"
ohdeerStyle["size"] = 1

deerStyle = NodeStyle()
deerStyle['hz_line_color'] = otherdeercolor
deerStyle['vt_line_color'] = singletoncolor
deerStyle['hz_line_width'] = widthsingl
deerStyle['vt_line_width'] = widthsingl
ohdeerStyle["fgcolor"] = "black"
ohdeerStyle["size"] = 1

generalStyle = NodeStyle()
generalStyle["fgcolor"] = "black"
generalStyle["size"] = 1

for leaves in tree:
	if "/deer/" in leaves.name:
		if clusterIDs in leaves.name:
			leaves.set_style(ohdeerStyle)
		else:
			leaves.set_style(deerStyle)
	else:
		leaves.set_style(generalStyle)

colorlist = ["#8dd3c7","#ffffb3","#bebada","#fb8072","#80b1d3","#fdb462","#b3de69","#fccde5","#d9d9d9","#bc80bd","#ccebc5","#ffed6f","#7fc97f","#beaed4","#fdc086","#ffff99","#386cb0","#f0027f", "#e08214","#fdb863","#fee0b6","#f7f7f7","#d8daeb","#b2abd2"]

j=0
for clid in clusters:
	#set styles for clusters
	stylename = "nsStyle"+str(clid)
	leafstylename ="lsStyle"+str(clid)
	exec(stylename+"=NodeStyle()")
	transp = "66"
	color = colorlist[j]+ transp #Additional digits are for transparency
	#print("%i %s" %(clid,color))
	exec(stylename+"['bgcolor'] = '"+color+"'")
	exec(stylename+"['hz_line_color'] = 'black'")
	exec(stylename+"['vt_line_color'] = 'black'")
	exec(stylename+"['vt_line_width'] = "+str(widthcluster))
	exec(stylename+"['hz_line_width'] = "+str(widthcluster))
	
	exec(leafstylename+"=NodeStyle()")
	exec(leafstylename+"['vt_line_width'] = "+str(widthcluster))
	exec(leafstylename+"['hz_line_width'] = "+str(widthcluster))
	exec(leafstylename+"['hz_line_color'] = 'black'")
	exec(leafstylename+"['vt_line_color'] = 'black'")
	
	###use it to set styles for clusters
	for n in tree.get_common_ancestor(clusters[clid]).traverse():
		exec("n.set_style(lsStyle"+str(clid)+")")
	commonnode = tree.get_common_ancestor(clusters[clid])
	exec("commonnode.set_style(nsStyle"+str(clid)+")")
	j+=1





#tree.show(tree_style=ts)

tree.render(outfile, w = 1200, units= 'px', dpi = 350, tree_style = ts)
