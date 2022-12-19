#!/usr/bin/env python

import os
import sys
import random
import ete3
import optparse
from Bio import AlignIO

parser=optparse.OptionParser()
parser.add_option('-t', '--tree', help='Treefile in newick format', type='str')
parser.add_option('-i', '--iterations', default = 3, help='number of iterations to increase the tree', type='int')
parser.add_option('-a', '--alignment', help='Input alignment in fasta format', type='str')
parser.add_option('-s', '--start', default = 0, help='gene start coordinate', type='int')
parser.add_option('-e','--end', default = 0, help='gene end coordinate', type ='int')
parser.add_option('-f','--informat',default = 0, help='tree format for input tree: default format is 0 format in ete3',type = 'int')
parser.add_option('-r','--outformat',default = 1, help='tree format for output tree: default format is 1 format in ete3',type = 'int')
parser.add_option('-o', '--outfolder', help='path to folder with trees and alignments', default = '', type='str')
parser.add_option('-n', '--nodes', help='clusters file', type='str')
parser.add_option('-c', '--cluster', help='cluster id', type='str')

subtreesize =120
minhumanleaves = 5

##get options
options, args=parser.parse_args()
tree = ete3.Tree(options.tree, format=options.informat) #input tree is from treetime and nonbinary
#levels = options.levels
outfolder = options.outfolder
clusterfile = options.nodes
start = options.start
end = options.end
iterations = options.iterations

ClusterDict = dict()
#Read clusters
with open(clusterfile, 'r') as clin:
	for l in clin:
		l=l.rstrip('\n')
		cluster, leaf = l.split('\t')
		if not cluster in ClusterDict:
			ClusterDict[cluster] = list()
		ClusterDict[cluster].append(leaf)

###############################
#alignment 
alignment = AlignIO.read(options.alignment, "fasta")

if (start > 0) & (end > 0):
	alignmentslice = alignment[:,start:end+1]
else:
	alignmentslice = alignment
	
#create folder
if not os.path.exists(outfolder):
    os.makedirs(outfolder)
else:
	print("WARNING: folder already exists")
   
cl = options.cluster
#get subtree
clustertree = tree.get_common_ancestor(ClusterDict[cl])
clustertree.name = '#1' #this line to create forground branch (I've placed it incorrectly before)

print(cl)
print(clustertree)
###
def check_num_humans(t):
	humcheck = -1
	humlist = list()
	for upleaves in t.get_leaves():
		if not "deer" in upleaves.name:
			humlist.append(upleaves.name)
			if len(humlist) > minhumanleaves:
				humcheck = 1
				break
	return(humcheck)
###
def do_check(t,c):
	if c == 1:
		twoup = t
	else:
		twoup = t.up
	return(twoup)
###
levelup = clustertree.up

for j in range(0,iterations):
	check = check_num_humans(levelup)
	if check == 1:
		break
	levelup = do_check(levelup, check)
	
#print(levelup)
	
#check if subsetting is required
allleaveslevel = levelup.get_leaves()
if len(allleaveslevel) > subtreesize:
	print("Large tree "+str(len(allleaveslevel)))
	random.shuffle(allleaveslevel)
	nondeerlist = list()
	for elem in allleaveslevel:
		if not "deer" in elem.name:
			nondeerlist.append(elem)
	numdelete = len(allleaveslevel) - subtreesize
	leavestodelete = nondeerlist[0:numdelete]
	for ld in leavestodelete:
		ld.delete()
	leavestokeep = list(set(allleaveslevel) - set(leavestodelete))
	allleaveslevel = leavestokeep
	print(len(levelup.get_leaves()))
#get names of leaves to keep
leavestokeepnames = list()
for lk in allleaveslevel:
	leavestokeepnames.append(lk.name)
########################
##subset alignment
renamedict = dict()
with open(outfolder+"/"+cl+'.fasta', 'w') as ouf:
	counter = 1
	for record in alignmentslice:
		newname = str(record.id).split("|")[0]
		parts = newname.split("/")
		phyname = ""
		if len(parts) > 1:
			if parts[1] == "deer":
				if parts[3].startswith("OH"):
					phyname = "dOH_"+str(counter)
				else:
					phyname = "d_"+str(counter)
			else:
				phyname = "h_"+str(counter)			
		else:
			phyname = parts[0][:7]
		renamedict[record.id] = phyname
		if record.id in leavestokeepnames:
			ouf.write(">%s\n%s\n" % (phyname, record.seq))
		counter += 1
		
#AlignIO.convert(outfolder+'/tmp.fasta', 'fasta', outfolder+"/"+cl+'.phy', 'phylip')
#clean up
#os.system('rm '+outfolder+'/tmp.fasta')		
##############################
#rename leaves on final tree and store names
with open(outfolder+"/"+cl+".leaves", 'w') as leavesf:
	finaltree = levelup
	for fl in finaltree.get_leaves():
		print(fl.name+"\t"+renamedict[fl.name])
		leavesf.write(fl.name+"\t"+renamedict[fl.name]+"\n")
		fl.name = renamedict[fl.name]
	finaltree.write(format=options.outformat, outfile=outfolder+"/"+cl+".nwk")

###get ctl files for codeml analysis
#those are the ctls for branch analysis 

with open(outfolder+"/"+cl+"_branch.ctl", 'w') as brctl, open(outfolder+"/"+cl+"_branch.fixed.ctl", 'w') as fixedctl:
	
	brctl.write("seqfile = "+outfolder+"/"+cl+'.fasta\n')
	fixedctl.write("seqfile = "+outfolder+"/"+cl+'.fasta\n')
	brctl.write("outfile = "+outfolder+"/"+cl+"_branch_amb.txt\n")
	fixedctl.write("outfile = "+outfolder+"/"+cl+"_branch_amb.fixed.txt\n")
	brctl.write("treefile = "+outfolder+"/"+cl+".nwk\n")
	fixedctl.write("treefile = "+outfolder+"/"+cl+".nwk\n")
	brctl.write("""noisy = 9
verbose = 1
runmode = 0
seqtype = 1
CodonFreq = 2
model = 2
NSsites = 2
Mgene = 0
icode = 0
fix_kappa = 0
kappa = 2
fix_omega = 0
omega = 1
RateAncestor = 0   * (0,1,2): rates (alpha>0) or ancestral states (1 or 2)
Small_Diff = .45e-6  * Default value.
cleandata = 0
""")
	fixedctl.write("""noisy = 9
verbose = 1
runmode = 0
seqtype = 1
CodonFreq = 2
model = 2
NSsites = 2
Mgene = 0
icode = 0
fix_kappa = 0
kappa = 2
fix_omega = 1
omega = 1
RateAncestor = 0   * (0,1,2): rates (alpha>0) or ancestral states (1 or 2)
Small_Diff = .45e-6  * Default value.
cleandata = 0
""")
