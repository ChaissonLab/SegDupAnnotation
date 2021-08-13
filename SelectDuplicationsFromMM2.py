#!/usr/bin/env python
import sys
alnFile=open(sys.argv[1])
alnType=sys.argv[2]
alns=[line.split() for line in alnFile]

i=0
d=0
isSam=False

if alnType == "paf":
    geneIdx=0
    endIdx=8
    startIdx=7
else:
    geneIdx=3
    endIdx=2
    startIdx=1

while i < len(alns):
    if alns[0][0] == "#":
        isSam=True
        i+=1
        continue
    j=i
#    if i == 27492:
#    import pdb
#    pdb.set_trace()
    while j < len(alns) and alns[i][geneIdx] == alns[j][geneIdx]:
        j+=1
    if j > i+1:
        # Found potential duplication

        if alnType == "paf":
            nMatchOrig=int(alns[i][10])
            origGeneLength=int(alns[i][1])
            ident=int(alns[i][9])/origGeneLength
        else:
            geneName=alns[i][3]
            origInterval=geneName.split("/")[-1]
            origStart=int(origInterval.split(":")[-1].split("-")[0])
            origEnd=int(origInterval.split(":")[-1].split("-")[1])
            origGeneLength=origEnd-origStart
            ident=float(alns[i][8])
        matchIdx=[i]
#        if "Fbxw14" in alns[i][0]:
#            import pdb
#            pdb.set_trace()
        dupAlnLength=int(alns[i][endIdx]) - int(alns[i][startIdx])
        alns[i].append(dupAlnLength/origGeneLength)
        alns[i].append(ident)
        for k in range(i+1,j):
            dupAlnLength=int(alns[k][endIdx]) - int(alns[k][startIdx])
            alns[k].append(str(dupAlnLength / origGeneLength))            
            if alnType == "paf":
                dupNMatch=int(alns[k][10])
                ident=int(alns[i][9])/origGeneLength                
#                for f in alns[k][12:]:
#                    if "dv:f" in f:
#                        ident=1-float(f.split(":")[-1])
                alns[k].append(ident)
                if dupAlnLength / origGeneLength > 0.9 and origGeneLength / dupAlnLength > 0.9 and ident > 0.9:
                    matchIdx.append(k)
            else:
                ident=float(alns[k][8])
                alns[k].append(ident)
                if dupAlnLength / origGeneLength > 0.9 and origGeneLength / dupAlnLength > 0.9 and ident > 0.9:
                    matchIdx.append(k)
                    
        if len(matchIdx) > 1:
            for k in matchIdx:                
                sys.stdout.write("\t".join(alns[k][0:12] + [str(v) for v in alns[k][-2:]]) + "\n")
        d+=1
    i=j
                
            
