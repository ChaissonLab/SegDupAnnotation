#!/usr/bin/env python
import sys
import pysam
import argparse
ap = argparse.ArgumentParser(description="Annotate repeat content of a bed file")
ap.add_argument("bed", help="Count fraction repeat of these regions")
ap.add_argument("ref", help="Regions are on this fasta file")

args=ap.parse_args()
bed=open(args.bed)
ref=pysam.FastaFile(args.ref)
for line in bed:
    vals=line.split()
    seq=ref.fetch(vals[0], int(vals[1]), int(vals[2]))
    n = seq.count("N") + seq.count("a") + seq.count("c") + seq.count("g") + seq.count("t") + seq.count("n")
    print("{:2.2f}".format(n/len(seq)))
