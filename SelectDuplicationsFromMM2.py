#!/usr/bin/env python
import sys
alnFile=open(sys.argv[1])
alns=[line.split() for line in alnFile]

i=0
d=0
while i < len(alns):
    j=i
#    if i == 27492:
#        import pdb
#        pdb.set_trace()
    while j < len(alns) and alns[i][0] == alns[j][0]:
        j+=1
    if j > i+1:
        # Found potential duplication
        origGeneLength=int(alns[i][8]) - int(alns[i][7])
        nMatchOrig=int(alns[i][10])
        matchIdx=[i]
#        if "Fbxw14" in alns[i][0]:
#            import pdb
#            pdb.set_trace()
        for k in range(i+1,j):
            dupAlnLength=int(alns[k][8]) - int(alns[k][7])
            dupNMatch=int(alns[k][10])
            if dupAlnLength / origGeneLength > 0.9 and origGeneLength / dupAlnLength > 0.9 and dupNMatch / dupAlnLength > 0.9:
                matchIdx.append(k)
        if len(matchIdx) > 1:
            for k in matchIdx:            
                sys.stdout.write("\t".join(alns[k]) + "\n")
        d+=1
    i=j
                
            
