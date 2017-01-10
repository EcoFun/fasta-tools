#!/usr/bin/env python2
# script to compute local GC content with sliding windows

__author__ = "Ludovic Duvaux"
__maintainer__ = "Ludovic Duvaux"
__license__ = "GPL_v3"

import sys, argparse
from Bio import SeqIO
from Bio.SeqUtils import GC

# 0.1) get options from commands lines
parser = argparse.ArgumentParser(description='Compute GC local in sliding windows')
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

if step > size:
	print "WARNING: window sliding step bigger than window size"

# if no sliding window
if step <1: step = size

#~print size, step
#~sys.exit()

#~finp = "../00_data/VeinB04.pacbio.fasta"
print "chromosome\tstart\tend\tPc_GC"
for rec in SeqIO.parse(finp, "fasta"):
	sta=0
	seq = rec.seq
	l = len(seq)
	while sta < l:
		sto = sta + size
		if sto > l: sto = l
		gc = GC(seq[sta:sto])
		v = [rec.name, sta+1, sto, gc]
		res = "\t".join([str(x) for x in v])
		print res
		# at the end increment sta by step
		sta += step
