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

    while j < len(vals) and Overlaps(vals[i], vals[j]):
        j+=1
    sys.stdout.write("\t".join(vals[i]) + "\n")
    i=j
