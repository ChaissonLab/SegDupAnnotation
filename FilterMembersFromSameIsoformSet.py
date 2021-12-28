#!/usr/bin/env python

import sys
if sys.argv[1] == "stdin":
    inFile=sys.stdin
else:
    inFile = open(sys.argv[1])

allVals = [line.split() for line in inFile.readlines()]

isoforms={}
geneIdx=None

def ParseGeneNameFromIDList(nameList):
    entries=nameList.split("|")
    if len(entries) >= 9:
        geneIdx=j
        accIdx=4
        isoformName=allVals[i][j].split("|")[4]
        geneName=allVals[i][j].split("|")[5]
        return (True, isoformName, geneName)
    elif len(entries) == 2:
        isoformName=nameList
        geneName=entries[0]
        return( True, isoformName, geneName)
    else:
        return (False, None, None)
    
        

for i in range(0,len(allVals)):
    #
    # First determine the index of the gene name for rules that are not bed12
    #
    accIdx=None
    
    if geneIdx is None:
        for j in range(0,len(allVals[i])):
            (foundGene, isoformName, geneName ) = ParseGeneNameFromIDList(allVals[i][j])
            if foundGene:
                geneIdx=j
                break
    if geneIdx is None:
        sys.stderr.write("ERROR, gene not found in " + "\t".join(allVals[i]) + "\n")
        sys.exit(1)
               
    (foundGene, isoformName, geneName) = ParseGeneNameFromIDList(allVals[i][geneIdx])
#    geneName=geneAcc[accIdx]
    if geneName not in isoforms:
        isoforms[geneName] = {}
    isoforms[geneName][isoformName] = True
    

firstIsoform = {}
for gene in isoforms.keys():
    firstIsoform[gene] = list(isoforms[gene].keys())[0]

for i in range(0,len(allVals)):
    
    
    (foundGene, isoformName, geneName) = ParseGeneNameFromIDList(allVals[i][geneIdx])
    if isoformName == firstIsoform[geneName]:
        sys.stdout.write("\t".join(allVals[i]) + "\n")

