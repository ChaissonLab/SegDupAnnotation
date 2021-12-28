#!/usr/bin/env python
import sys
genes={}
inFile=open(sys.argv[1])

header=inFile.readline()
vals=header.split()
hmmCountIndex=None
geneIndex=0
resolvedIndex=None
for i in range(0,len(vals)):

    if vals[i] == "copy":
        hmmCountIndex=i
    if vals[i] == "gene":
        geneIndex=i
    if vals[i] == "resolved":
        resolvedIndex=i
if hmmCountIndex is None:
    print("ERROR could not find column with copy")
    sys.exit(1)

print("Gene\tasmCount\tcolCount")
for line in inFile:
    vals=line.split()
    gene=vals[geneIndex]
    if gene not in genes:
        genes[gene] = [0,0]
        
    if vals[resolvedIndex] == "multi":
        genes[gene][0]+=1
    else:
        genes[gene][1]+=int(vals[hmmCountIndex])
for gene in genes:
    print(gene + "\t" + str(genes[gene][0]) + "\t" + str(genes[gene][1]))
    
