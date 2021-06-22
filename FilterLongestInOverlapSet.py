#!/usr/bin/env python
import sys
if sys.argv[1] == "stdin":
    inFile=sys.stdin
else:
    inFile = open(sys.argv[1])


allVals = [line.split() for line in inFile.readlines()]

isoforms={}
i=0

def GetForwardCoordinates(aln):
#    pdb.set_trace()
    return (int(aln[1]), int(aln[2]))
        
def FractionOverlap(alnsA, alnsB):
    if alnsA[0] != alnsB[0]:
        return 0
    
    fa=GetForwardCoordinates(alnsA)
    fb=GetForwardCoordinates(alnsB)
    ovp=0
    la=fa[1]-fa[0]
    lb=fb[1]-fb[0]
    if (fa[0] <= fb[0] and fa[1] >= fb[0]):
        ovp=min(fb[1], fa[1]) - fb[0]
    if (fa[0] >= fb[0] and fa[0] <= fb[1]):
        ovp = min(fb[1], fa[1]) - fa[0]

    return max(ovp/la, ovp/lb)



while i < len(allVals):
    j=i
    while j < len(allVals):
        ovp= FractionOverlap(allVals[i], allVals[j])
        if ovp < 0.5:
            break
        j+=1
    maxLen=0
    maxK=0
    for k in range(i,j):
        l=int(allVals[i][2]) - int(allVals[i][1])
        if l > maxLen:
            maxLen=l
            maxK=k
    sys.stdout.write("\t".join(allVals[k]) + "\n")
    i=j
