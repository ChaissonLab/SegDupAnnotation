#!/usr/bin/env python
import sys
genes={}
header=sys.stdin.readline()
vals=header.split()
hmmCountIndex=None

for i in range(0,len(vals)):
    if vals[i] == "copy":
        hmmCountIndex=i
        break

if hmmCountIndex is None:
    print("ERROR could not find column with copy")
    sys.exit(1)

print("Gene\tasmCount\tcolCount")
for line in sys.stdin:
    vals=line.split()
    gene=vals[3]

    if gene not in genes:
        genes[gene] = [0,0]
    if vals[6] == "multi":
        genes[gene][0]+=1
    else:
        genes[gene][1]+=int(vals[hmmCountIndex])
for gene in genes:
    print(gene + "\t" + str(genes[gene][0]) + "\t" + str(genes[gene][1]))
    
