#!/usr/bin/env python
import sys
inFile = open(sys.argv[1])
minScore = int(sys.argv[2])
lines= [v.split() for v in inFile.readlines()]

for i in range(0,5):
    sys.stdout.write("\t".join(lines[i]) + "\n")
i=5
while i < len(lines):
    j=i

    while j < len(lines) and len(lines[i]) >9 and len(lines[j]) > 9 and lines[i][9] == lines[j][9] :
        j+=1
    topScore = int(lines[i][0])
    if topScore > minScore:
        k=i+1
        sys.stdout.write("\t".join(lines[i]) + "\n")
        while k < j and int(lines[k][0]) / topScore > 0.8:
            sys.stdout.write("\t".join(lines[k]) + "\n")
            k+=1
    i=j

    
