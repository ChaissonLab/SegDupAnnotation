#!/usr/bin/env python
import sys
import argparse 

ap = argparse.ArgumentParser(description="Use a graph based approach for counting duplications")

ap.add_argument("bed", help="Input bed file, for sensitivity, use paired-approach.")
ap.add_argument("condensed", help="Condensed graph")
ap.add_argument("tab", help="Duplication table")
ap.add_argument("--overlap", help="Fraction overlap", type=float, default=0.5)
ap.add_argument("--full", help="Full output graph", default=None)

args     = ap.parse_args()

inFile   = open(args.bed)
outFile  = args.condensed
outTableName = args.tab
vals = [l.split() for l in inFile]
keys=[]
import networkx as nx

def FractionOverlap(alnsA, alnsB):
    if alnsA[0] != alnsB[0]:
        return 0
    
    ovp=0
    la=alnsA[2]-alnsA[1]
    lb=alnsB[2]-alnsB[1]
    
    if (alnsA[1] <= alnsB[1] and alnsA[2] >= alnsB[1]):
        ovp=min(alnsB[2], alnsA[2]) - alnsB[1]
    if (alnsA[1] >= alnsB[1] and alnsA[1] <= alnsB[2]):
        ovp = min(alnsB[2], alnsA[2]) - alnsA[1]
    l=max(la,lb)

    ovp=ovp/l
    return ovp


class Node:
    def __init__(self, key):
        self.key=key
        self.adj = {}
        self.traversed = False

    def HasEdge(self, q):
        return q in self.adj


groups = {}
i=0
g=0
index=0
groupIndex=[-1]*len(vals)
vertices = {}
black=0
red=1

#
# G is the final graph, and bG is the modified (contracted graph)
G=nx.Graph()
bG=nx.Graph()
left=0
right=1
allKeys={}
for i in range(0,len(vals)):
    if vals[i][0][0] == "#":
        continue
    key="_".join(vals[i][0:3])
    keyl="_".join(vals[i][0:3])+"_"+str(i) #+"_L"

    G.add_node(keyl, chrom=vals[i][0], start=int(vals[i][1]), end=int(vals[i][2]), side=left, index=i)
    bG.add_node(keyl, chrom=vals[i][0], start=int(vals[i][1]), end=int(vals[i][2]), side=left, index=i)
    keyr="_".join(vals[i][3:6])+"_" + str(i) # + "_R"
    G.add_node(keyr, chrom=vals[i][0], start=int(vals[i][1]), end=int(vals[i][2]), side=left, index=i)
    bG.add_node(keyr, chrom=vals[i][0], start=int(vals[i][1]), end=int(vals[i][2]), side=left, index=i)
    allKeys[keyl] = True



    
for i in range(0,len(vals)):
    if vals[i][0][0] == "#":
        continue
    key="_".join(vals[i][0:3]) + "_" + str(i) # + "_L"
    destKey="_".join(vals[i][3:6]) + "_" + str(i) #+ "_R"

    G.add_edge(key, destKey, color=red)
       

    if i not in groups:
        groups[i] = g
        g+=1
    j=i+1
    if j >= len(vals):
        continue
    chromi=vals[i][0]
    starti=int(vals[i][1])
    endi=int(vals[i][2])
    leni=endi-starti;
    chromj=vals[j][0]
    startj=int(vals[j][1])
    endj=int(vals[j][2])

    while (j < len(vals) and chromi == chromj and startj < endi ):
#        if starti==11998287 or startj == 11998287:
#            import pdb
#            pdb.set_trace()

        lenj=endj-startj;

        overlap=FractionOverlap((chromi, starti, endi), (chromj, startj, endj))
        destKey = "_".join(vals[j][0:3])+ "_" + str(j)# + "_L"
        
        if (overlap > args.overlap):
            G.add_edge(key, destKey, color=black)
            bG.add_edge(key, destKey, color=black)

#            AddAdj(vertices, key, destKey, black)
            groups[j] = groups[i]
        j+=1
        if (j < len(vals)):
            chromj=vals[j][0]
            startj=int(vals[j][1])
            endj=int(vals[j][2])
        index +=1
        if (index % 50000 == 0):
            sys.stderr.write(str(i) + "/" + str(len(vals)) + "\t" + str(j) + "\t" + str(g) + "\n")
    i+=1


if (args.full is not None):
    nx.write_gml(G, outFile)

# Now contract nodes.
nRemoved=0
nAdded=0
curIndex=0
compIndex=[-1]*len(vals)
compSize=[1]*len(vals)

for comp in nx.connected_components(G):
    compList=list(comp)
    s=len(compList)
    for node in compList:
        idx=G.nodes[node]["index"]
        compIndex[idx] = curIndex
        compSize[idx] = s
    curIndex+=1

outTable=open(outTableName, 'w')
for i in range(0,len(compSize)):
    outTable.write(str(compSize[i]) + "\t" + str(compIndex[i]) + "\n")

for n in G.nodes:
    G.nodes[n]["index"] = str(G.nodes[n]["index"])

nx.write_gml(G, args.condensed)
compIndex=0


#outTable.write("\n".join([str(i) for i in compSize]))


