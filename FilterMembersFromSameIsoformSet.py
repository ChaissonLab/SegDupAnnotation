#!/usr/bin/env python

import sys
if sys.argv[1] == "stdin":
    inFile=sys.stdin
else:
    inFile = open(sys.argv[1])

allVals = [line.split() for line in inFile.readlines()]

isoforms={}
geneIdx=None
for i in range(0,len(allVals)):
    #
    # First determine the index of the gene name for rules that are not bed12
    #
    if geneIdx is None:
        for j in range(0,len(allVals[i])):
            if len(allVals[i][j].split("|")) >= 9:
                geneIdx=j
                break
    if geneIdx is None:
        sys.stderr.write("ERROR, gene not found in " + "\t".join(allVals[i]) + "\n")
        sys.exit(1)
               
    geneAcc=allVals[i][geneIdx].split("|")
    geneName=geneAcc[5]
    if geneName not in isoforms:
        isoforms[geneName] = {}
    isoforms[geneName][geneAcc[4]] = True
    

firstIsoform = {}
for gene in isoforms.keys():
    firstIsoform[gene] = list(isoforms[gene].keys())[0]

for i in range(0,len(allVals)):    
    geneAcc=allVals[i][geneIdx].split("|")
    isoform=geneAcc[4]
    geneName=geneAcc[5]
    if isoform == firstIsoform[geneName]:
        sys.stdout.write("\t".join(allVals[i]) + "\n")

