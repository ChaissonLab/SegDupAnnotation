#!/usr/bin/env python
import sys

lines = [ l.split() for l in sys.stdin ]

genes = {}
excluded = 0
for line in lines:
    gencode=line[3].split("|")
    gene=gencode[5]
    if gencode[5] in genes:
        excluded+=1
        continue
    else:
        genes[gene] = True
        sys.stdout.write("\t".join(line) + "\n")

