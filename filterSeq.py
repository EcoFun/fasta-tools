#!/usr/bin/python
from Bio import SeqIO
import sys

"""
This script remove sequences from a fasta file according to a list of id
"""

def main(argv = None):
  if argv is None:
    argv = sys.argv
  if len(argv) <= 2:
    print "Usage: filterSeq.py <fasta_file> <id_file>"
    return 1

  fastaFileName = argv[1]
  idFileName = argv[2]
  #print "Fasta file:", fastaFileName
  #print "Id file:", idFileName

  # Read id into list
  ids = set()
  for l in open(idFileName):
    ids.add(l.strip())
  for seq_record in SeqIO.parse(fastaFileName, "fasta"):
    if seq_record.id in ids:
      pass
    else:
      print seq_record.format("fasta"),
  return 0

if __name__ == "__main__":
  sys.exit(main())
