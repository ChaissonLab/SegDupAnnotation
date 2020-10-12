#!/usr/bin/env python

import sys
import argparse
ap=argparse.ArgumentParser(description="Filter pacbio alignments to extract alignments from longest mapped subread.")
ap.add_argument("--input", help="input file.", default="/dev/stdin")
args=ap.parse_args()

reads = {}
inFile = open(args.input)
nRead=0
for line in inFile:
    vals=line.split()
    subread=vals[3]
    readBase="/".join(vals[3].split("/")[0:2])
    if readBase not in reads:
        reads[readBase] = {}
    if subread not in reads[readBase]:
        reads[readBase][subread] = []
    reads[readBase][subread].append(vals)
    nRead+=1
    if nRead % 100000 == 0:
        sys.stderr.write("read " + str(nRead/1000) + "k alignments\n")





def GetForwardCoordinates(aln):
#    pdb.set_trace()
    return (int(aln[5]), int(aln[6]))
#    if aln[4] == "0":
#
#    if aln[4] == "1":
#        s=[int(i) for i in aln[3].split("/")[2].split("_")]
#        l = s[1] - s[0]
#        start = l-int(aln[5])
#        end=l-int(aln[6])+1
#        return (start,end)
#
def FractionOverlap(alnsA, alnsB):
    fa=GetForwardCoordinates(alnsA)
    fb=GetForwardCoordinates(alnsB)
    ovp=0
    la=fa[1]-fa[0]
    lb=fb[1]-fb[0]
    if (fa[0] <= fb[0] and fa[1] >= fb[0]):
        ovp=min(fb[1], fa[1]) - fb[0]
    if (fa[0] >= fb[0] and fa[0] <= fb[1]):
        ovp = min(fb[1], fa[1]) - fa[0]

    return max(ovp/la, ovp/lb)

    
for readBase in reads.keys():
    # Compute length of each alignment
    subreads=reads[readBase]
    alnLengths = []
    iSets = []
    for subread in subreads:
        alns = subreads[subread]
        alnLengths=[ int(alns[i][2]) - int(alns[i][1]) for i in range(0,len(alns)) ]
        index=[v[1] for v in sorted(list(zip(alnLengths, range(0,len(alns)))))]
        iSet = [index[0]]
        for i in index[1:]:
            foundOverlap=False
            for j in range(0,len(iSet)):
                if FractionOverlap(alns[i], alns[iSet[j]]) > 0.5:
                    foundOverlap=True
            if foundOverlap is False:
                iSet.append(i)

        alnLength=0
        for  i in iSet:
            alnLength+=int(alns[i][2]) - int(alns[i][1])
        alnLengths.append(alnLength)
        iSets.append(iSet)

    keyLengths=sorted(list(zip(alnLengths, subreads.keys(), iSets)), reverse=True)

    subread=keyLengths[0][1]
    iset=keyLengths[0][2]
                     
    for i in iset:
        sys.stdout.write("\t".join(subreads[subread][i]) + "\n")
    
        

        
    
