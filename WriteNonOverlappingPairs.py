#!/usr/bin/env python
import sys
vals = [ l.strip().split() for l in sys.stdin ]
i = 0
curGene = 0
def Overlaps(a,b, astart=0, bstart=0):
    aStart= int(a[1+astart])
    aEnd  = int(a[2+astart])
    aChrom = a[0+astart]
    bStart= int(b[1+bstart])
    bEnd  = int(b[2+bstart])
    bChrom = b[0+bstart]
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

    #
    # Collect all lines that overlap on the first interval
    while j < len(vals) and Overlaps(vals[i], vals[j]):
        j+=1
        
    exclude=[False]* (j-i)
    for ii in range(i,j-1):
        for jj in range(ii+1,j):
            if (Overlaps(vals[ii], vals[jj], astart=3,bstart=3)):
                exclude[ii-i] = True
    for ii in range(i,j):
        if (exclude[ii-i] == False and Overlaps(vals[ii], vals[ii], bstart=3) == False):
            sys.stdout.write("\t".join(vals[ii]) + "\n")
    i=j
