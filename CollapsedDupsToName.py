#!/usr/bin/env python


import sys
inFile=open(sys.argv[1])
outFile=open(sys.argv[2], 'w')

for line in inFile:
    vals=line.split()
    vals[3] = vals[3].split("|")[5]
    outFile.write("\t".join(vals) + "\n")

