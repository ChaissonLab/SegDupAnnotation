#!/usr/bin/env python

import sys
import argparse
ap=argparse.ArgumentParser(description="Filter pacbio alignments to extract alignments from longest mapped subread.")
ap.add_argument("--input", help="input file.", default="/dev/stdin")




def coordinate(r_len, orient, start, end):
    if orient == 0:
        f=(start, end, end-start)
    elif orient == 1:
        r_s = r_len - end
        r_e = r_len - start
        f=(r_s, r_e, r_e - r_s)
    else:
        raise ValueError('Found orient: {}'.format(orient))
    return f


def getstart(elem):
    l1 = elem
    rb = l1[3].split("/")[2].split("_")
    r_len = int(rb[1])-int(rb[0])
    c1 = coordinate( r_len, int(l1[4]),int(l1[5]),int(l1[6]) )
    return c1[0]


def span(alignmentList, non):
    sorted_intervals = []
    span = 0
    for i in range(0, len(alignmentList) ) :
        if i in non:
            sorted_intervals.append((0, 0))
        else:
            l1 = alignmentList[i]
            rb = l1[3].split("/")[2].split("_")
            r_len = int(rb[1])-int(rb[0])
            c1 = coordinate(r_len, int(l1[4]), int(l1[5]), int(l1[6]))
            #print("i"+str(i))
            sorted_intervals.append((c1[0], c1[1]))
    low = sorted_intervals[0][0]
    high = sorted_intervals[0][1]
    for v in range(1, len(sorted_intervals) ):
        if sorted_intervals[v][0] <= high:  # new interval overlaps current run
            #print(sorted_intervals[v][1], high)
            if v in non:
                continue
            if high > sorted_intervals[v][1]:
                continue
            else:
                high = sorted_intervals[v][1] # merge with the current run
        else:  # current run is over
            span += (high-low)
            low, high = sorted_intervals[v][0], sorted_intervals[v][1]  # start new run
    span += (high-low)
    return span  # end the final run

def NonOverlappingAlignedBases(alignmentList):
    p=len(alignmentList)
    non=set()
    cs=()
    for r in range(0,p-1):
        l1= alignmentList[r]
        rb= l1[3].split("/")[2].split("_")
        r_len = int(rb[1])-int(rb[0])
        c1=coordinate(r_len, int(l1[4]),int(l1[5]),int(l1[6]))
        #span=span+c1[2]
        l2= alignmentList[r+1]
        c2=coordinate(r_len, int(l2[4]),int(l2[5]),int(l2[6]))
        if c1[0]==c2[0] and c1[1]==c2[1]:
                non.add(i)
                non.add(j)
                continue

    cs = ( span(alignmentList, non), non)
    return cs



args=ap.parse_args()
inFile = open(args.input)
reads={}
#
#read in all alignments and their subreads
#
for line in inFile:
    line=line.rstrip()+"\n"
    vals=line.split("\t")

    if (int(vals[7]) < 30):
        continue

    readName=vals[3]
    readVals=readName.split("/")
    molecule="/".join(readVals[0:2])
    if molecule not in reads:
        reads[molecule] = {}
    subread = "/".join(readVals[0:3])
    if subread not in reads[molecule]:
        reads[molecule][subread] = []
    reads[molecule][subread].append(vals)


for molecule in reads.keys():
    subreadKeys = reads[molecule].keys()
    y= len(reads[molecule])
    if y == 1:
        for key in subreadKeys:
            u=len(reads[molecule][key])
            for i in range(0,u):
                sys.stdout.write("\t".join(reads[molecule][key][i]))
        continue
    c=()
    max=0
    index=None
    for subreadAlignmentSet in subreadKeys:
        z=len(reads[molecule][subreadAlignmentSet])
        ls=reads[molecule][subreadAlignmentSet]
        s_ls=sorted(ls,key=getstart)
        c = NonOverlappingAlignedBases( s_ls )
        compare=int(c[0])
        #print("comp"+str(compare))
        if compare > max :
            index=subreadAlignmentSet
            max=compare
    #print("max"+str(max))
    for j in range(0,len(reads[molecule][index])  ):
        #print(c[1])
        if j not in c[1] :
            sys.stdout.write("\t".join(reads[molecule][index][j]) )

