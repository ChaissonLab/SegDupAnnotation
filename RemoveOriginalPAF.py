#!/usr/bin/env python
import sys
import re

def Overlap(s1, e1, s2, e2):
    if e1 < s2 or s1 > e2:
        return 0
    
    maxS=max(s1,s2)
    minE=min(e1,e2)
    return (minE-maxS)/(e2-s2)
 
annotate=False
if sys.argv[1] == "stdin" or sys.argv[1] == "-" or sys.argv[1] == "/sys/stdin":
    inFile = sys.stdin
else:
    inFile = open(sys.argv[1])
    if len(sys.argv) > 2:
        if sys.argv[2] == "annotate":
            annotate=True
for line in inFile:
    vals=line.split()
    chrom=vals[5]
    start=int(vals[7])
    end=int(vals[8])
    nameVals=vals[0].split("/")[-1]
    nameRgn=re.split('(.*):(.*)-(.*)', nameVals)
    srcChrom=nameRgn[1]
    srcStart=int(nameRgn[2])
    srcEnd=int(nameRgn[3])
    ovp=Overlap(start, end, srcStart, srcEnd)

    if (chrom == srcChrom and ovp > 0.5):
        status="Original"
        if annotate:            
            sys.stdout.write(line.strip() + "\tOriginal\n")
    else:
        if annotate:
            sys.stdout.write(line.strip() + "\tCopy\n")
        else:
            sys.stdout.write(line)
