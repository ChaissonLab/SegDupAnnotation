#!/usr/bin/env python
import sys
import argparse
import pysam

ap=argparse.ArgumentParser(description="Remove sedef output if one alignment has too much masked")
ap.add_argument("--genome", help="Genome", required=True)
ap.add_argument("--pct", help="Min percent masked", required=True, type=float)
ap.add_argument("--sedef", help="Sedef output", required=True)
ap.add_argument("--maskedBed", help="Write a bed file of all the masked entries here.", required=True)
args=ap.parse_args()

sedef = open(args.sedef)
asm=pysam.FastaFile(args.genome)
mb=open(args.maskedBed, "w")

for line in sedef:
    vals=line.split()
    rgn1=vals[0] + ":" + vals[1] + "-" + vals[2]
    rgn2=vals[3] + ":" + vals[4] + "-" + vals[5]
    seq1=asm.fetch(vals[0], int(vals[1]), int(vals[2]))
    seq2=asm.fetch(vals[3], int(vals[4]), int(vals[5]))
    mask1=seq1.count("a") + seq1.count("c") + seq1.count("g") + seq1.count("t")
    mask2=seq2.count("a") + seq2.count("c") + seq2.count("g") + seq2.count("t")
    if mask1/len(seq1) > args.pct or mask2/len(seq2) > args.pct:
        mb.write(vals[0] + "\t" + vals[1] + "\t" + vals[2] + "\n")
        mb.write(vals[3] + "\t" + vals[4] + "\t" + vals[5] + "\n")
        continue
    sys.stdout.write(line)
    
