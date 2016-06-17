#!/bin/env python2
# usage:
    # ./compare_cleaning_methods.py raw_outliers_files columns_header sample_file

__author__ = "Ludovic Duvaux"
__maintainer__ = "Ludovic Duvaux"
__license__ = "GPL_v3"

import sys, numpy, argparse
v=0

# 0.1) get options from commands lines
parser = argparse.ArgumentParser(description='Sort alphabetically sequences from a fasta file')
parser.add_argument('fasta_file', help='Files of the raw lists of outliers', nargs=1)
parser.add_argument('-f', help='Name of a sequence to be displayed first', nargs=1)

# 0.2) set options
args = parser.parse_args()
#~print args

s = args.fasta_file
#~print s
fsam = s[0]
#~print fsam

# 1) set initial dictionary
res = {}
with open(fsam) as sam:
    for l in sam:
        #~print l
        if l[0] ==">":
            key = l.strip("> ").strip()
            res[key] = [l.strip()]
        else:
            res[key].append(l.strip())
#~print res

# print final ali
    # print first sequence
if args.f != None:
    first = args.f[0]
    keys = res.keys()
    m = [s for s in keys if first in s ]
    if len(m) == 0:
        print("ERROR: Incorrect string for option -f")
        exit ("       No corresponding match\n")
    elif  len(m) > 1:
        print("ERROR: Ambiguous string for option -f")
        exit ("       Too many possibilities\n")
    
    fir = m[0]
    ll = res[fir]
    for s in ll:
        print s
    del res[fir]    # remove that entry from list

# print the rest
sorted_key = sorted(res.keys(), key=lambda s: s.lower())
#~print sorted_key
for i in sorted_key:
    ll = res[i]
    for s in ll:
        print s
