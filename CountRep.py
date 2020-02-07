#!/usr/bin/env python
import pysam
import argparse 

ap = argparse.ArgumentParser(description="Count how much a target sequence is repetitive")

ap.add_argument("--bed", help="Input bed file", required=True)
ap.add_argument("--ref", help="Reference", required=True)

args=ap.parse_args()

dups=open(args.bed)
ref=pysam.FastaFile(args.ref)
dups.readline()
for line in dups:
    vals=line.split()
    sc=vals[0]
    ss=int(vals[1])
    se=int(vals[2])
    tc=vals[3]
    ts=int(vals[4])
    te=int(vals[5])

    acc=float(vals[7])
    
    if (acc > 40):
        continue
    if (abs((se-ss)  - (te-ts)) < 200):
        sSeq=ref.fetch(vals[0], ss, se)
        tSeq=ref.fetch(vals[3], ts, te)
        nls=sSeq.count("a") + sSeq.count("c") + sSeq.count("g") + sSeq.count("t")
        nlt=tSeq.count("a") + tSeq.count("c") + tSeq.count("g") + tSeq.count("t")

        fs = nls/(se-ss)
        ft = nlt/(te-ts)
        print("{}\t{}\t{}\t{:2.2f}\t{}\t{}\t{}\t{:2.2f}".format(sc, ss, se, fs*100, tc, ts, te, ft*100))

    
