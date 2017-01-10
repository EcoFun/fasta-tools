#!/bin/env python2

__author__ = "Ludovic Duvaux"
__maintainer__ = "Ludovic Duvaux"
__license__ = "GPL_v3"

### I) packages and functions
import argparse, collections, numpy, sys
from os.path import basename
from Bio import SeqIO

### II) get options from commands lines
# II.1) define options
parser = argparse.ArgumentParser(description='### Count the character frequency per sequence from an alignment')
parser.add_argument('align', help='Alignment file', nargs=1)

# II.2) parse options
args = parser.parse_args()
fin = args.align[0]

# III) SCRIPT
# 1) open file and count characters
sys.stderr.write("Compute character frequency (may take a while)\n")
res = {}
handle = open(fin, "rU")
for record in SeqIO.parse(handle, "fasta"):
    key = record.id
    ct = collections.Counter(record.seq.upper())
    res[key] = ct

handle.close()

# 2) get all possible values of characters
sorted_key = sorted(res.keys(), key=lambda s: s.lower())
good = ['A', 'C', 'G', 'T']
all_val = [] + good
for k in sorted_key:
    v = res[k].keys()
    for i in v:
        if i not in all_val:
            all_val += i


# 3) write result file
sys.stderr.write("### Write output\n")
print "Sample\t" + "\t".join([ str(s) for s in all_val ]) + "\tNo_ACGT"
for i in sorted_key:
    ct = res[i]
    # sum different ACGT
    som = 0
    for v in all_val:
        if v not in good:
            som += int(ct[v])
    
    vec = i + "\t" + "\t".join([ str(ct[s]) for s in all_val ])
    vec += "\t" + str(som)
    print vec
