#!/usr/bin/env python
import sys
inFile=open(sys.argv[1])
for line in inFile:
    if len(line) > 0 and line[0] == "#":
        continue
    vals=line.split()
    if vals[6] != "PASS":
        continue
    chrom=vals[0]
    start=vals[1]
    kvp=vals[7].split(";")
    end=kvp[1].split("=")[1]
    cn=vals[-1].split(":")[0]
    sys.stdout.write(chrom+"\t" + start + "\t" + end + "\t" + cn + "\n")
