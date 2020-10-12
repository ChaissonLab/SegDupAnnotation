#!/usr/bin/env python
import pysam
import sys
inFile = pysam.AlignmentFile(sys.argv[1])
outFile = open(sys.argv[2],'w')
multiMaps={}
for gene in inFile:
    name = gene.query_name
    tup=gene.cigartuples
    if tup is None:
        continue
    clip=0
    i=0

    while (i < len(tup) and (tup[i][0] == 4  or tup[i][0] == 5)):
        clip+= tup[i][1]
        i+=1
    j=len(tup)-1
    while (j > i and (tup[j][0] == 4  or tup[j][0] == 5)):
        clip+= tup[j][1]
        j-=1

    if gene.query_sequence is None:
        continue
    transcriptLen=len(gene.query_sequence)
    if float(clip)/transcriptLen < 0.2:
        print(gene.to_string())
        if name not in multiMaps:
            multiMaps[name] = []

        multiMaps[name].append(gene)



for geneName in multiMaps:
    if len(multiMaps[geneName]) > 1:
        i=0
        for g in multiMaps[geneName]:
            outFile.write(str(i) + "\t" + g.to_string() + "\n")
            i+=1
