#!/bin/env python2

__author__ = "Ludovic Duvaux"
__maintainer__ = "Ludovic Duvaux"
__license__ = "GPL_v3"

### I) packages and functions
import sys, numpy, argparse
from os.path import basename

def entries_from_files (fil):
    lis = []
    if fil != None:
        fil = fil[0]
        with open(fil) as f:
            for l in f:
                lis.append(l.strip())
    return lis;

def update_table (list_files, res, bad, col, og, ingp=True, keys=[]):
    if (ingp is not True) and (ingp is not False):
        print "ERROR: ingp must be False or True"
        return
    for fil in list_files:
        with open(fil) as f:
            # get dic key
            fil = basename(fil)
            sa = fil.split('.')
            key = [s for s in sa if "BUSCO" in s][0]

            for l in f:
                if l[0] == ">":
                    # detect indiv
                    l = l.strip().strip("> ")
                    ind = l.split("_")[0]
                    # update dic
                    if ingp:
                        # update only if 'ind' not in outg nor bad guys 
                        if ind not in (bad + og):
                            res[key][col] += 1
                    elif not ingp:
                        # update only if 'ind' in outgroups
                        if ind in og:   
                            res[key][col] += 1
                    else:
                        print "INTERNAL ERROR"
                        return
    return res;


### II) get options from commands lines
# II.1) define options
parser = argparse.ArgumentParser(description='Count the number of OD-seq outliers per locus.')
parser.add_argument('-e', '--exclude', help='File of individuals to be ignored in the analysis (i.e. low quality genomes that will be discarded from downstream analyses).', nargs=1)    # 0.1.2) path to samples to be excluded
parser.add_argument('-f', '--file-outgp', help='file specifying outgroups if OD-seq analyses haven\'t been performed separately for the ingroup and the outgroups (not compatible with the \'-o\' flag).', nargs=1)
parser.add_argument('-i', '--ingp_fils', help='Results files of OD-seq analyses. If used with the \'-o\' flag, specify results file for the ingroup only.', nargs='+', required=True)
parser.add_argument('-o', '--outgp_fils', help='Results files of OD-seq analyses for the outgroups only (not compatible with the \'-f\' flag).', nargs='+')
parser.add_argument('-v', '--verbose', help='outgroups separated by commas (useful if only one bunch of OD-seq results to analyse).', nargs=1, default=-1)   # verbosity

# II.2) parse options
args = parser.parse_args()
v = args.verbose
if v != -1: v = int(v[0])
fin = args.ingp_fils
fout = args.outgp_fils
foutg = args.file_outgp
fexc = args.exclude

if v == 0: print args;
if v == 1: print v;
if v == 2: print fin
if v == 3: print fout
if v == 4: print foutg

if fout != None and foutg != None:
    exit ("ERROR: the -f and -o flags cannot be used together")
if fout == None: fout = fin

# II.3) get file entries
bad = entries_from_files(fexc)  # samples to be ignored
og = entries_from_files(foutg)  # list outgroups from file
if v == 5: print bad
if v == 5: print og


### III) SCRIPT
# 1) set initial dictionary of loci -> count outliers
res = {}
for f in fin:
    # get dic key
    f = basename(f)
    sa = f.split('.')
    gn = [s for s in sa if "BUSCO" in s][0]
    # set dic
    res[gn] = [0, 0, 0]
sorted_key = sorted(res.keys(), key=lambda s: s.lower())
if v == 6: print res

# 3.1) count number of outlier locus from ingroup
res = update_table (fin, res, bad, col=0, og=og, ingp=True, keys=sorted_key)
if v == 7: print res

# 3.2) count number of outlier locus from outgroup
res = update_table (fout, res, bad, col=1, og=og, ingp=False, keys=sorted_key)
if v == 8: print res

# 4) sum
for i in res.keys():
    res[i][2] = sum(res[i][:2])
if v == 9: print res

# 5) print output
print "Locus\tNingroup\tNoutgroup\tTotal"

for i in sorted_key:
    val = "\t".join([ str(s) for s in res[i] ])
    print "%s\t%s" % (i, val)

