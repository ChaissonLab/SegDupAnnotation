#!/usr/bin/env python
import sys
inFile = open(sys.argv[1])
vals = [l.split() for l in inFile]

groups = {}
i=0
g=0
index=0
while i < len(vals):
    if i not in groups:
        groups[i] = g
        g+=1
    j=i
    starti=int(vals[i][1])
    endi=int(vals[i][2])
    leni=endi-starti;
    startj=int(vals[j][1])
    endj=int(vals[j][2])
    
    while (j < len(vals) and startj < endi and startj < starti+leni/4):
        lenj=endj-startj;
        overlap=(endi-startj)/leni
        if (overlap > 0.8):
            groups[j] = groups[i]
        j+=1

        startj=int(vals[j][1])
        endj=int(vals[j][2])
        index +=1
        if (index % 50000 == 0):
            sys.stderr.write(str(i) + "/" + str(len(vals)) + "\t" + str(j) + "\t" + str(g) + "\n")
    i+=1

groupCount = {}
repetitive = {}
for i in groups:
    if groups[i] not in groupCount:
        groupCount[i] = 0
    groupCount[i] += 1
    if groupCount[i] > 20:
        repetitive[i] = True

for i in range(0,len(vals)):
    g=groups[i]
    if repetitive[i]:
        sys.stdout.write("\t".join(vals[i]) + "\n")



        
    
