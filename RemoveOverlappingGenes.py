#!/usr/bin/env python
import sys
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


prevGene=None
inFile=sys.stdin
lines = inFile.readlines()

i=0
j=0
keep=[True] * len(lines)
while i < len(lines):
    vals=lines[i].split()
    curGene=vals[3]
    j=i+1
    while (j < len(lines)):
        nextVals=lines[j].split()
        nextGene=nextVals[3]
        if nextGene == curGene:
            if FractionOverlap(vals, nextVals) > 0.1:
                keep[j] = False
            j+=1
        else:
            break
    i+=1

for i in range(0,len(lines)):
    if keep[i]:
        sys.stdout.write(lines[i])


