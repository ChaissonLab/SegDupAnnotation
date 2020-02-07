#!/usr/bin/env python

import sys
inFile = open(sys.argv[1])

allVals = [line.split() for line in inFile.readlines()]

isoforms={}
for i in range(0,len(allVals)):
    geneAcc=allVals[i][3].split("|")
    geneName=geneAcc[5]
    if geneName not in isoforms:
        isoforms[geneName] = {}
    isoforms[geneName][geneAcc[4]] = True
    

firstIsoform = {}
for gene in isoforms.keys():
    firstIsoform[gene] = list(isoforms[gene].keys())[0]

for i in range(0,len(allVals)):
    geneAcc=allVals[i][3].split("|")
    isoform=geneAcc[4]
    geneName=geneAcc[5]
    if isoform == firstIsoform[geneName]:
        sys.stdout.write("\t".join(allVals[i]) + "\n")

