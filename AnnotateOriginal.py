#!/usr/bin/env python
import sys
import re
inFile = open(sys.argv[1])
for line in inFile:
    vals=line.split()
    chrom=vals[0]
    start=int(vals[1])
    end=int(vals[2])
    nameVals=vals[3].split("/")[-1]
    nameRgn=re.split('(.*):(.*)-(.*)', nameVals)
    srcChrom=nameRgn[1]
    srcStart=int(nameRgn[2])
    srcEnd=int(nameRgn[3])
    
    if (chrom == srcChrom and (abs(srcStart - start) < 100 and abs(srcEnd-end) < 100) or ( ( start >= srcStart and  start < srcEnd) or (end >= srcStart and end < srcEnd) or (start <= srcStart and end >= srcEnd) ) ):
        sys.stdout.write(line.strip() + "\tOriginal\n")
    else:
        sys.stdout.write(line.strip() + "\tCopy\n")
