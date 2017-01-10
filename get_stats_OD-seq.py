#!/bin/env python2

__author__ = "Ludovic Duvaux"
__maintainer__ = "Ludovic Duvaux"
__license__ = "GPL_v3"

import sys, numpy
v="false"

# collect data
ls =[]
for l in sys.stdin:
    l=l.strip("> ").strip("\n")
    l=l.split("BUSCO")[0]
    if v == "true": print l
    ls.append(l)

if v == "true": print(ls)

# get stats
ct=numpy.unique(ls, return_counts=True)
print len(ct[0])

# print results 
for i in range(0,len(ct[0])-1):
    print "%s\t%s" % (ct[0][i], ct[1][i])
