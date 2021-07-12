#!/usr/bin/env python
import sys
mapFileName = sys.argv[1]
if mapFileName == "NO_OP":
    for line in sys.stdin:
        sys.stdout.write(line)
else:
    mapFile=open(mapFileName)
    nameMap={}
    for line in mapFile:
        vals=line.split()
        nameMap[vals[0]] = vals[1]
    for line in sys.stdin:
        vals=line[0]
        for i in range(0,len(vals)):
            if vals[i] in nameMap:
                vals[i] = nameMap[vals[i]]
            elif "var" in vals[i] and vals[i][3:] in nameMap:
                vals[i] = "var" + nameMap[vals[i][3:]]
        sys.stdout.write("\t".join(vals) + "\n")
        
