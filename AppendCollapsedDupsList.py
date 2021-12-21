#!/usr/bin/env python
import sys
alnDupCnFile=open(sys.argv[1])
depthDupCnFile=open(sys.argv[2])
header=alnDupCnFile.readline()

aln = [ l.split() for l in alnDupCnFile ]
depth = [ l.split() for l in depthDupCnFile ]

alnGenes = { a[3] : True  for a in aln }
sys.stdout.write(header.rstrip() + "\tresolved\n")
for i in range(0, len(aln)):
    sys.stdout.write("\t".join(aln[i] + ["multi"]) + "\n")
    
for d in depth:
    extra=str(max(0,int(d[-1])-2))
    d[-1] = "1"
    d.append("ORIG")
    d.append(extra)
    d.append("collapse")
    sys.stdout.write("\t".join(d) + "\n")



