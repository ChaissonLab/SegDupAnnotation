#!/usr/bin/env python
import sys
alnFile=open(sys.argv[1])
alns=[line.split() for line in alnFile]

i=0
d=0
isSam=False

geneIdx=3
endIdx=2
startIdx=1
copyNumberIdx=13

while i < len(alns):
    if alns[0][0] == "#":
        isSam=True
        i+=1
        continue
    j=i
    while j < len(alns) and alns[i][geneIdx] == alns[j][geneIdx]:
        j+=1
    if j > i+1:
        # Found potential duplication

        geneName=alns[i][3]
        origInterval=geneName.split("/")[-1]
        origStart=int(origInterval.split(":")[-1].split("-")[0])
        origEnd=int(origInterval.split(":")[-1].split("-")[1])
        origGeneLength=origEnd-origStart
        ident=float(alns[i][8])
            
    matchIdx=[i]
    if j > i+1:
        #
        # There are resolved copies of a gene.  Add the identity here.
        #
        dupAlnLength=int(alns[i][endIdx]) - int(alns[i][startIdx])
        alns[i].append(dupAlnLength/origGeneLength)
        alns[i].append(ident)
        for k in range(i+1,j):
            dupAlnLength=int(alns[k][endIdx]) - int(alns[k][startIdx])
            alns[k].append(str(dupAlnLength / origGeneLength))            
            ident=float(alns[k][8])
            alns[k].append(ident)
#            import pdb
#            pdb.set_trace()
#            print("{}\t{}\t{:2.2f}\t{:2.2f}\t{:2.2f}".format(i,j,dupAlnLength / origGeneLength, origGeneLength / dupAlnLength, ident))
            if dupAlnLength / origGeneLength > 0.9 and origGeneLength / dupAlnLength > 0.9 and ident > 0.9:
                matchIdx.append(k)
    else:
        #
        # Collapsed duplication, add 1's for 100% identity.
        #
        alns[i].append(1)
        alns[i].append(1)
    if len(matchIdx) > 1 or int(float(alns[i][copyNumberIdx]) ) > 2:
        for k in matchIdx:
            sys.stdout.write("\t".join(alns[k][0:12] + [alns[k][13]] + [str(v) for v in alns[k][-2:]]) + "\n")
        d+=1
    i=j
                
            
