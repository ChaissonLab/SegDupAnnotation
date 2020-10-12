#!/usr/bin/env python
import sys
trfFile=open(sys.argv[1])
trfBed=open(sys.argv[2],'w')
contig=None
for line in trfFile:
    if len(line) > 0 and line[0] == "@":
        contig=line.rstrip()[1:]
    else:
        vals=line.split()
        trfBed.write(contig + "\t" + vals[0] + "\t" + vals[1] + "\n")
