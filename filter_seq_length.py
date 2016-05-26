#!/bin/env python2

__author__ = "Ludovic Duvaux"
__maintainer__ = "Ludovic Duvaux"
__license__ = "GPL_v3"

from sys import argv, exit

usage="""
SYNOPSIS:
'{0}' fasta_file length output_file

DESCRIPTION:
Filter sequences from a fasta file (e.g. contigs) based on length.

ARGUMENTS:
    fasta_file          input fasta file
    length              threshold below which sequences are discarded
    output_file         name of the output fasta file

WARNING: the script work only with fasta files having sequences on only 
one line!
""".format(basename(argv[0]))

if len(argv) != 4:
    print "ERROR: wrong number of arguments (there should be 3)"
    exit(usage)

inp = argv[1]       # inpout
leng = int(argv[2]) # length filter
out = argv[3]       # output


o = open(out, "w")
print "# Write file %s" % out

i = 0
with open(inp) as f:
    for l in f:
        if l[0] != '>':
            print "Sequence on multiple lines!"
            exit(usage)
        else
            m = next(inp)   # record following line (should be sequence line)
            lg = len(m)
            if lg > leng:
                o.write(l + "\n")
                o.write(m + "\n")
            else
                print "skip sequence %s (length=%s)" % lg

o.write('\n')
f.close()
o.close()
