#!/usr/bin/env python2

__author__ = "Ludovic Duvaux"
__maintainer__ = "Ludovic Duvaux"
__license__ = "GPL_v3"

from Bio import SeqIO
import sys, argparse

# 0.1) get options from commands lines
parser = argparse.ArgumentParser(description='Format unusual fasta styles into mainstream fasta. If a list of individual IDs is provided, filter out those individuals.')
parser.add_argument('fasta_file', help='Multifasta alignement to be formatted', nargs=1)
parser.add_argument('-i', '--id_file', help='List file of IDs to be removed', nargs=1)

argv = parser.parse_args()
fastaFileName = argv.fasta_file[0]
#~sys.stderr.write(fastaFileName + "\n")

# Read id into list
if argv.id_file:
	idFileName = argv.id_file[0]
	#~sys.stderr.write(idFileName + "\n")
	ids = []
	with open(idFileName) as f:
		for l in f:
			ids.append(l.strip())
	#~sys.stderr.write("\t".join(ids) + "\n")

# write
with open(fastaFileName) as f:
	for seq_record in SeqIO.parse(f, "fasta"):
		sys.stderr.write("Process " + seq_record.id + "\n")
		if argv.id_file and seq_record.id in ids:
			#~sys.stderr.write("PASS %s\n" % seq_record.id)
			continue
		print seq_record.upper().format("fasta").strip()
