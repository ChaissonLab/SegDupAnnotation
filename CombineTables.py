#!/usr/bin/env python
import sys

nArgs=len(sys.argv)
nSamp=int((nArgs-1)/2)


files = [open(fn) for fn in sys.argv[1:(nSamp+1)]]
names = sys.argv[(1+nSamp):]

for f in files:
    f.readline()
tables = [ { l.split()[0] : l.split()[1] for l in f } for f in files ]

allGenes={}
for tab in tables:
    for k in tab.keys():
        allGenes[k] = ['1']*nSamp
print("gene\t" + "\t".join(names))
for gene in allGenes.keys():
    for i in range(0,len(tables)):
        if gene in tables[i]:
            allGenes[gene][i] = tables[i][gene]
    print(gene + "\t" + "\t".join(allGenes[gene]))


