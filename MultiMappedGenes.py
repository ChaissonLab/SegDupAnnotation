#!/usr/bin/env python
import sys
allPairwise = []
for line in sys.stdin:
    vals=line.split()
    allPairwise.append(vals)

i=0
toFilter = {}
toKeep = {}
while i < len(allPairwise):
    geneStart = i
    geneEnd = i+1
    overlappedGenes={}
    #
    # Check how many different genes overlap current gene.
    allGenes={}
    curGene = allPairwise[geneStart][3]
    while geneEnd < len(allPairwise) and allPairwise[geneEnd][3] == allPairwise[geneStart][3]:
        allGenes[allPairwise[geneEnd][15]] = True
        allGenes[allPairwise[geneEnd][3]]  = True
        if allPairwise[geneEnd][15] != allPairwise[geneEnd][3]:
            overlappedGenes[allPairwise[geneEnd][15]] = True
        geneEnd+=1
    i=geneEnd
    if curGene not in toFilter:
        toKeep[curGene] = True
    if len(overlappedGenes) > 0:
        for g in overlappedGenes.keys():
            if g not in toKeep:
                toFilter[g] = True

for f in toFilter.keys():
    print(f)


        

        
        
        
        
    
