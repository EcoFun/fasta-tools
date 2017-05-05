#!/usr/bin/env python2.7
# script to compute local GC content with sliding windows

__author__ = "Ludovic Duvaux"
__maintainer__ = "Ludovic Duvaux"
__license__ = "GPL_v3"

import sys, argparse
from Bio import SeqIO
from Bio.SeqUtils import GC

# 0.1) get options from commands lines
parser = argparse.ArgumentParser(description='Compute GC local in sliding windows')
parser.add_argument('-b', '--bed', 
	help='bed file specifying windows coordinates for which to compute GC content (BEWARE: bed have 0-based coordinates)',
	nargs=1, default=[None])
parser.add_argument('-s', '--size', 
	help='Windows size in which to compute the GC in bp [1000]',
	nargs=1, default=[1000], type=int)
parser.add_argument('-S', '--step', 
	help='Sliding step of the windows in bp. By default the sliding step is "size" [size]',
	nargs=1, default=[0], type=int)
parser.add_argument('fasta', help='multi fasta file', nargs=1)

# 0.2) set options
args = parser.parse_args()
#~print args
size = int(args.size[0])
step = int(args.step[0])
finp = args.fasta[0]
bed = args.bed[0]

if step > size:
	print "WARNING: window sliding step bigger than window size"

# 1) compute GC if bed provided
if bed:
	# load bed file in dict
	#~i=0
	targets = {}
	with open(bed) as b:
		for l in b:
			if l[0]=="#":
				continue
			ll = l.strip().split()	# bed file, so 0-based coord like python
			try:
				targets[ll[0]].append([int(ll[1]), int(ll[2])]) 
			except KeyError:
				targets[ll[0]] = [ [int(ll[1]), int(ll[2])] ]
		
	# then start reading fasta file
	print "chromosome\tstart\tend\tPc_GC"
	for rec in SeqIO.parse(finp, "fasta"):
		seq = rec.seq
		scaff = rec.name
		if scaff not in targets.keys():
			continue
		else:
			for coord in targets[scaff]:
				sta = coord[0]	# bed file, so 0-based coord like python
				sto = coord[1]	# bed file, so 0-based coord like python
				gc = GC(seq[sta:(sto+1)])	# bed file, so 0-based coord like python
				v = [rec.name, sta+1, sto+1, gc]	# bed file, so 0-based coord like python
				res = "\t".join([str(x) for x in v])
				print res
				#~i += 1

# 2) else compute GC for the whole genome
else:
	if step <1: 
		step = size	# if no sliding window
	print "chromosome\tstart\tend\tPc_GC"
	for rec in SeqIO.parse(finp, "fasta"):
		sta=0	# beware the coordinates!
		seq = rec.seq
		l = len(seq)
		while sta < l:
			sto = sta + size
			if sto > l: sto = l
			gc = GC(seq[sta:sto])	# beware the coordinates!
			v = [rec.name, sta+1, sto, gc]
			res = "\t".join([str(x) for x in v])
			print res
			# at the end increment sta by step
			sta += step


