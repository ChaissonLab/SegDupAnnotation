#!/usr/bin/env python
import sys

vals = [ l.strip().split() for l in sys.stdin ]
i = 0
curGene = 0
def Overlaps(a,b):
    aStart= int(a[1])
    aEnd  = int(a[2])
    aChrom = a[0]
    bStart= int(b[1])
    bEnd  = int(b[2])
    bChrom = b[0]
    if aChrom != bChrom:
        return False
    if aStart >= bStart  and aStart <= bEnd:
        return True
    if aEnd >= bStart and aEnd <= bEnd:
        return True
    if aStart <= bStart and aEnd >= bEnd:
        return True
    return False
    
while i < len(vals):
    j=i

    while i < len(vals) and j < len(vals) and vals[j][3] == vals[i][3]:
        j+=1
    leni= int(vals[i][2]) - int(vals[i][1])
    ovpStr = "\t".join(vals[i]) + "\n"
    nOvp = 0
    for k in range(i+1,j):
        lenk = int(vals[k][2]) - int(vals[k][1])
        if lenk / leni > 0.8:
            noOverlap = True
            for k2 in range(i,k):
                if Overlaps(vals[k2], vals[k]) == True:
                    noOverlap = False
            if noOverlap:
                ovpStr += "\t".join(vals[k]) + "\n"
                nOvp+=1

    if nOvp > 0:
        sys.stdout.write(ovpStr)

    i=j
        
