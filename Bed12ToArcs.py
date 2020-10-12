#!/usr/bin/env python
import sys
dupGenes={}
arcsFile=open(sys.argv[1],'w')
namesFile=open(sys.argv[2],'w')
for line in sys.stdin:
    vals=line.split()
    if vals[3] not in dupGenes:
        dupGenes[vals[3]] = []
    dupGenes[vals[3]].append(vals)


for dupGene in dupGenes.keys():
    l=len(dupGenes[dupGene])
    if l > 1:
        for i in range(0,l-1):
            for j in range(i+1,l):
                vi=dupGenes[dupGene][i]
                vj=dupGenes[dupGene][j]

                nExon=max(int(vi[9]), int(vj[9]))
                if vi[0] == vj[0]:
                    if nExon == 1:
                        col="paired-6-qual-1"
                    else:
                        col="paired-6-qual-2"
                else:
                    if nExon == 1:
                        col="paired-6-qual-3"
                    else:
                        col="paired-6-qual-4"
                arcsFile.write("var" + vi[0] + "\t" + vi[1] + "\t" + vi[2] + "\tvar" + vj[0] + "\t" + vj[1] + "\t" + vj[2]+ "\tcolor="+col + "\n")
        for i in range(0,l):
            vi=dupGenes[dupGene][i]
            namesFile.write("var" + vi[0] + "\t" + vi[1] + "\t" + vi[2] + "\t" + dupGene + "\n")


