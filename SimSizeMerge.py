#!/usr/bin/env python
import sys
prevLine=[]
nMerged=0
cIntv = None
import pdb
inFile = open(sys.argv[1])
for line in inFile:
    if line[0] == "#":
        continue
    vals=line.split()
    if float(vals[7]) > 20:
        continue
    intv=[vals[0], int(vals[1]), int(vals[2])]
    isOvp=False
    if cIntv is not None and cIntv[2] == 3337039 and intv[2] == 3336633:
        print(intv)
        pdb.set_trace()

    if cIntv is not None:
        if vals[0] == cIntv[0] and \
           ((cIntv[1]<= intv[1] and cIntv[2]> intv[1]) or \
            (cIntv[2] > intv[1] and cIntv[2] <= intv[2]) or \
            (cIntv[1] <= intv[1] and cIntv[2] > intv[2])):
            #
            # The two intervals overlap
            #
            lc = cIntv[2]- cIntv[1]
            li = intv[2] - intv[1]
            if (abs(lc-li) < 0.2*max(lc,li)):
                cIntv[1] = min(cIntv[1], intv[1])
                cIntv[2] = max(cIntv[2], intv[2])
                isOvp=True
                nMerged +=1
        if isOvp == False:
            if cIntv is not None and cIntv[2] == 3337039:
                print(intv)
                pdb.set_trace()

            print(cIntv[0] + "\t" + str(cIntv[1]) + "\t" + str(cIntv[2]) + "\t" + str(nMerged))
            cIntv = intv
            nMerged = 1
    else:
        cIntv = intv
        
        
