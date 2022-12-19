#!/usr/bin/env python

import sys
import optparse

parser=optparse.OptionParser()
parser.add_option('-i', '--infile', help='', type='str')
parser.add_option('-o', '--outfile', help='', type='str')

options, args=parser.parse_args()

Dict =dict()

with open(options.infile, 'r') as inf:
	for l in inf:
		l=l.rstrip('\n')
		p=l.rsplit(' ', 1)
		if not p[0] in Dict:
			Dict[p[0]]=list()
		Dict[p[0]].append(p[1])
		
#print(Dict)


with open(options.outfile, 'w') as ouf:
	for k in Dict:
		ouf.write('%s %i %s\n' %(k, len(Dict[k]), ",".join(Dict[k])))
