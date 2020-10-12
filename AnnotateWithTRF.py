#!/usr/bin/env python

import argparse
import sys
import intervaltree as itree

ap = argparse.ArgumentParser(description="Annotate a gap bed file with an associated TRF table.")
ap.add_argument("trf", help="Input TRF table file.")
ap.add_argument("bed_out", help="Output BED file.")
args = ap.parse_args()

trfFile = open(args.trf, 'r')
prevSeqName = ""

annotations = {}
for line in trfFile:
    if (line[0] == '@'):
        seqName = line[1:].strip()
    else:
        vals = line.split()
        if (seqName not in annotations):
            annotations[seqName] = itree.IntervalTree()
        annotations[seqName].addi(int(vals[0]), int(vals[1]))




bedFile = open(args.bed, 'r')
bedOutFile = open(args.bed_out, 'w')
h = {}
for line in bedFile:
    if line[0] == "#":
        header = line.split()
        h={header[i] : i for i in range(0,len(header))}
        bedOutFile.write("nTR\tfracTR\n")
        continue
    vals = line.split()
    seqTitle = '/'.join(vals[0:3])
    totalTR = 0    
    if (seqTitle in annotations):
        seq = vals[h["svSeq"]]
        i = 0
        annotations[seqTitle].merge_overlaps()

        for intv in annotations[seqTitle]:
            totalTR += intv[1] - intv[0]
            seq = seq[:intv[0]] + seq[intv[0]:intv[1]].lower() + seq[intv[1]:]
        vals[5] = seq

        vals= [str(totalTR), "{:2.2f}".format(float(totalTR)/len(seq))]
    else:
        vals=["0","0"]

    bedOutFile.write('\t'.join(vals)  + "\n")
        
