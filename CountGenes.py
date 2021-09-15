#!/usr/bin/env python
import sys
genes={}
for line in sys.stdin:
    vals=line.split()
    gene=vals[3]

    if gene not in genes:
        genes[gene] = [0,0]
    if vals[6] == "multi":
        genes[gene][0]+=1
    else:
        genes[gene][1]+=int(vals[5])
for gene in genes:
    print(gene + "\t" + str(genes[gene][0]) + "\t" + str(genes[gene][1]))
    
