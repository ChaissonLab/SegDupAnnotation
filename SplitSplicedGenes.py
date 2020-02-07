#!/usr/bin/env python
import sys
if len(sys.argv) != 4:
    sys.stdout.write("Usage: SplitSplicedGenes.py inFile.bed multiExon.bed singleExon.bed\n")
    sys.exit(1)

inFile = open(sys.argv[1])
multiExon = open(sys.argv[2], 'w')
singleExon = open(sys.argv[3], 'w')
for line in inFile:
    vals=line.split()
    spliceSites=vals[10].split(",")
    sys.stderr.write(str(spliceSites) + "\n")
    if len(spliceSites) == 1:
        singleExon.write(line)
    else:
        multiExon.write(line)
