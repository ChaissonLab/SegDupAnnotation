#!/usr/bin/env python
import sys
dupGenes={}

arcsFile=open(sys.argv[1],'w')
namesFile=open(sys.argv[2],'w')
for line in sys.stdin:
    vals=line.split()
    if vals[0] not in dupGenes:
        dupGenes[vals[0]] = []
    dupGenes[vals[0]].append(vals)


for dupGene in dupGenes.keys():
    l=len(dupGenes[dupGene])
    if l > 1:
        for i in range(0,l):
            vi=dupGenes[dupGene][i]
            if vi[1] == vi[8]:
                col="paired-6-qual-1"
            else:
                col="paired-6-qual-3"                
            arcsFile.write("var" + vi[1] + "\t" + vi[2] + "\t" + vi[3] + "\tvar" + vi[8] + "\t" + vi[9] + "\t" + vi[10]+ "\tcolor="+col + "\n")
        for i in range(0,l):
            vi=dupGenes[dupGene][0]
            namesFile.write("var" + vi[1] + "\t" + vi[2] + "\t" + vi[3] + "\t" + dupGene + "\n")


