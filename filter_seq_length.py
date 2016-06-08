#!/bin/env python2

__author__ = "Ludovic Duvaux"
__maintainer__ = "Ludovic Duvaux"
__license__ = "GPL_v3"

from sys import argv, exit
from os.path import basename

usage="""
SYNOPSIS:
'{0}' fasta_file length output_file

DESCRIPTION:
Remove sequences shorter than 'length' from a fasta file.
Accept gz files.

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

inp = argv[1]       # input
leng = int(argv[2]) # length filter
out = argv[3]       # output

test_gz = inp[-2:]=="gz"
#~print test_gz
j = 1
if test_gz:
    import gzip
    if out[-2:]!="gz": out=out+".gz"
    o = gzip.open(out, "w")
    print "# Write file %s (compressed)" % out

    i = 0
    with gzip.open(inp) as f:
        for l in f:
            if l[0] != '>':
                print l
                print "Sequence on multiple lines!"
                exit(usage)
            else:
                m = next(f)   # record following line (should be sequence line)
                j += 2
                lg = len(m)
                if lg > leng:
                    o.write(l)
                    o.write(m)
                else:
                    print "skip sequence %s (length=%s)" % (l, lg)

else:
    o = open(out, "w")
    print "# Write file %s (uncompressed)" % out
    with open(inp) as f:
        for l in f:
            if l[0] != '>':
                print l
                print "Sequence on multiple lines!"
                exit(usage)
            else:
                m = next(f)   # record following line (should be sequence line)
                j += 2
                lg = len(m)
                if lg > leng:
                    o.write(l)
                    o.write(m)
                else:
                    print "skip sequence %s (length=%s)" % (l, lg)

print "File originally included %s sequences" % ((j-1)/2)
f.close()
o.close()
