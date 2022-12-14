#!/usr/bin/env python
import sys
alnFile=open(sys.argv[1])
alns=[line.split() for line in alnFile]
keep=[ True ] * len(alns)

endIdx=2
startIdx=1

for i in range(0,len(alns)):
    geneName=alns[i][3] # "A1CF|1/scaffold_10:35364605-35393550"
    origInterval=geneName.split("/")[-1] # "scaffold_10:35364605-35393550"
    origStart=int(origInterval.split(":")[-1].split("-")[0]) # "35364605"
    origEnd=int(origInterval.split(":")[-1].split("-")[1]) # "35393550"
    origGeneLength=origEnd-origStart # "28945"
    dupAlnLength=int(alns[i][endIdx]) - int(alns[i][startIdx]) # "28946"

    ident=float(alns[i][8]) # "report accuracy" from rule MappedSamIdentity (% match in bps of regions in assembly that map to gene model)
    if ident < 0.9 or origGeneLength / dupAlnLength < 0.9 or dupAlnLength / origGeneLength < 0.9:
        keep[i] = False
alns2=[]
alns3=[] # 2D array with removed stuff b/c low identity or length

for i in range(0,len(alns)):
    if keep[i]:
        alns2.append(alns[i])
    else:
        alns3.append(alns[i])
alns=alns2

i=0
d=0 # counts total duplications
isSam=False

geneIdx=3
copyNumberIdx=13

while i < len(alns):
    if alns[0][0] == "#":
        isSam=True
        i+=1
        continue
    j=i
    while j < len(alns) and alns[i][geneIdx] == alns[j][geneIdx]: # if gene names match (the whole "A1CF|1/scaffold_10:35364605-35393550")
        j+=1
    if j > i+1:
        # Found potential duplication

        geneName=alns[i][3]
        origInterval=geneName.split("/")[-1]
        origStart=int(origInterval.split(":")[-1].split("-")[0])
        origEnd=int(origInterval.split(":")[-1].split("-")[1])
        origGeneLength=origEnd-origStart
        ident=float(alns[i][8])
    
    matchIdx=[i] # list of indexes of matches, initially starts with itself, if resolved, additional indices added as appropriate
    noMatchIdx=[]
    noMatchIdxC=[]
    if j > i+1:
        #
        # There are resolved copies of a gene.  Add the identity here.
        #
        dupAlnLength=int(alns[i][endIdx]) - int(alns[i][startIdx])
        alns[i].append(dupAlnLength/origGeneLength) # col 14
        alns[i].append(ident) # col 15
        for k in range(i+1,j):
            dupAlnLength=int(alns[k][endIdx]) - int(alns[k][startIdx])
            alns[k].append(str(dupAlnLength / origGeneLength))  # col 14          
            ident=float(alns[k][8])
            alns[k].append(ident) # col 15
#            print("{}\t{}\t{:2.2f}\t{:2.2f}\t{:2.2f}".format(i,j,dupAlnLength / origGeneLength, origGeneLength / dupAlnLength, ident))
            if dupAlnLength / origGeneLength > 0.9 and origGeneLength / dupAlnLength > 0.9 and ident > 0.9:
                matchIdx.append(k)
            else:
                noMatchIdx.append(k)
    else:
        #
        # Collapsed duplication, add 1's for 100% identity.
        #
        alns[i].append(1) # col 14
        alns[i].append(1) # col 15
        noMatchIdxC.append(i)
    if len(matchIdx) > 1 or int(float(alns[i][copyNumberIdx]) ) > 2: # if dup found, print
        for k in matchIdx:
            sys.stdout.write("\t".join(alns[k][0:12] + [alns[k][13]] + [str(v) for v in alns[k][-2:]] + [alns[k][14]]) + "\n")
        d+=1
    else:
        for k in matchIdx:
            if alns[k][14] == ".":
                sys.stdout.write("\t".join(alns[k][0:12] + [alns[k][13]] + [str(v) for v in alns[k][-2:]] + ["filtOut:MappedSamIdentityDups-1"]) + "\n")
            else:
                sys.stdout.write("\t".join(alns[k][0:12] + [alns[k][13]] + [str(v) for v in alns[k][-2:]] + [alns[k][14]]) + "\n")
        for k in noMatchIdx:
            if alns[k][14] == ".":
                sys.stdout.write("\t".join(alns[k][0:12] + [alns[k][13]] + [str(v) for v in alns[k][-2:]] + ["filtOut:MappedSamIdentityDups-2"]) + "\n")
            else:
                sys.stdout.write("\t".join(alns[k][0:12] + [alns[k][13]] + [str(v) for v in alns[k][-2:]] + [alns[k][14]]) + "\n")
        #for k in noMatchIdxC: # Actually these collapses are already accounted for
        #    if alns[k][14] == ".":
        #        sys.stdout.write("\t".join(alns[k][0:12] + [alns[k][13]] + [str(v) for v in alns[k][-2:]] + ["filtOut:MappedSamIdentityDups-5"]) + "\n")
        #    else:
        #        sys.stdout.write("\t".join(alns[k][0:12] + [alns[k][13]] + [str(v) for v in alns[k][-2:]] + [alns[k][14]]) + "\n")
    i=j

i=0
while i < len(alns3):
    if alns3[0][0] == "#":
        isSam=True
        i+=1
        continue
    j=i
    while j < len(alns3) and alns3[i][geneIdx] == alns3[j][geneIdx]: # if gene names match (the whole "A1CF|1/scaffold_10:35364605-35393550")
        j+=1
    if j > i+1:
        # Found potential duplication

        geneName=alns3[i][3]
        origInterval=geneName.split("/")[-1]
        origStart=int(origInterval.split(":")[-1].split("-")[0])
        origEnd=int(origInterval.split(":")[-1].split("-")[1])
        origGeneLength=origEnd-origStart
        ident=float(alns3[i][8])
            
    matchIdx=[i] # list of indexes of matches, initially starts with itself, if resolved, additional indices added as appropriate
    noMatchIdx=[]
    noMatchIdxC=[]
    if j > i+1:
        #
        # There are resolved copies of a gene.  Add the identity here.
        #
        dupAlnLength=int(alns3[i][endIdx]) - int(alns3[i][startIdx])
        alns3[i].append(dupAlnLength/origGeneLength) # col 14
        alns3[i].append(ident) # col 15
        for k in range(i+1,j):
            dupAlnLength=int(alns3[k][endIdx]) - int(alns3[k][startIdx])
            alns3[k].append(str(dupAlnLength / origGeneLength))  # col 14          
            ident=float(alns3[k][8])
            alns3[k].append(ident) # col 15
            if dupAlnLength / origGeneLength > 0.9 and origGeneLength / dupAlnLength > 0.9 and ident > 0.9:
                matchIdx.append(k)
                print("aaa\t" + str(len(matchIdx)))
            else:
                noMatchIdx.append(k)
                #print("bbb\t" + str(len(noMatchIdx)))
    else:
        #
        # Collapsed duplication, add 1's for 100% identity.
        #
        alns3[i].append(1) # col 14
        alns3[i].append(1) # col 15
        noMatchIdxC.append(i)
    
    for k in matchIdx:
        if alns[k][14] == ".":
            sys.stdout.write("\t".join(alns3[k][0:12] + [alns3[k][13]] + [str(v) for v in alns3[k][-2:]] + ["filtOut:MappedSamIdentityDups-3"]) + "\n")
        else:
            sys.stdout.write("\t".join(alns3[k][0:12] + [alns3[k][13]] + [str(v) for v in alns3[k][-2:]] + [alns3[k][14]]) + "\n")
    for k in noMatchIdx:
        if alns[k][14] == ".":
            sys.stdout.write("\t".join(alns3[k][0:12] + [alns3[k][13]] + [str(v) for v in alns3[k][-2:]] + ["filtOut:MappedSamIdentityDups-4"]) + "\n")
        else:
            sys.stdout.write("\t".join(alns3[k][0:12] + [alns3[k][13]] + [str(v) for v in alns3[k][-2:]] + [alns3[k][14]]) + "\n")
    #for k in noMatchIdxC: # Actually these collapses are already accounted for
    #    if alns[k][14] == ".":
    #        sys.stdout.write("\t".join(alns3[k][0:12] + [alns3[k][13]] + [str(v) for v in alns3[k][-2:]] + ["filtOut:MappedSamIdentityDups-6"]) + "\n")
    #    else:
    #        sys.stdout.write("\t".join(alns3[k][0:12] + [alns3[k][13]] + [str(v) for v in alns3[k][-2:]] + [alns3[k][14]]) + "\n")

    i=j
