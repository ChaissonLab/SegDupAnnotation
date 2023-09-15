#!/usr/bin/env python
import sys

covFile = open(sys.argv[1])
bedFile = open(sys.argv[2])
avgCov  = float(sys.argv[3])
outFile = open(sys.argv[4], 'w')

genes = {}
for line in covFile:
    vals=line.split()
    genes[vals[0]] = vals[1]

for line in bedFile:
    vals=line.split()
    gene = vals[3]
    if gene in genes:
        if genes[gene] != "0":        
            outFile.write(line.strip() +  "\t{:2.2f}\t{}\n".format(float(genes[gene])/float(avgCov), genes[gene]))
        else:
            outFile.write(line.strip() + "\t0\t0\n")
    else:
        outFile.write(line.strip() + "\t0\n")
