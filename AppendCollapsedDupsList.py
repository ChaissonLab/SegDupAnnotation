#!/usr/bin/env python
import sys
alnDupCnFile=open(sys.argv[1])
depthDupCnFile=open(sys.argv[2])

aln = [ l.split() for l in alnDupCnFile ]
depth = [ l.split() for l in depthDupCnFile ]

alnGenes = { a[3] : True  for a in aln }
for i in range(0, len(aln)):
    sys.stdout.write("\t".join(aln[i] + ["multi"]) + "\n")
for d in depth:
    if d[3] not in alnGenes:
        extra=str(max(0,int(d[-1])-2))
        d[-1] = extra
        sys.stdout.write("\t".join(d + ["single"]) + "\n")



