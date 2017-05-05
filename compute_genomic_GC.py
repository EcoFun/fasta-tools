#!/usr/bin/env python2
# script to compute local GC content with sliding windows

__author__ = "Ludovic Duvaux"
__maintainer__ = "Ludovic Duvaux"
__license__ = "GPL_v3"

import sys, argparse
from Bio import SeqIO
from Bio.SeqUtils import GC

# 0.1) get options from commands lines
parser = argparse.ArgumentParser(description='Compute GC for a whole genome')
	nargs=1, default=[0], type=int)
parser.add_argument('fasta', help='multi fasta file', nargs=1)

# 0.2) set options
args = parser.parse_args()
finp = args.fasta[0]

for rec in SeqIO.parse(finp, "fasta"):
	seq = rec.seq
	
