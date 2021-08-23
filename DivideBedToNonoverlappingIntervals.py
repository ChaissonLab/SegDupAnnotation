#!/usr/bin/env python
import sys
curEnd=None
curChrom=None


from intervaltree import Interval, IntervalTree

trees={}
lineNumber=0
inFile=open(sys.argv[1])
for line in inFile:
    vals=line.split()
    if vals[0] not in trees:
        trees[vals[0]] = IntervalTree()
    lineNumber+=1
    trees[vals[0]].addi(int(vals[1]), int(vals[2]), lineNumber)



for chrom in trees.keys():
    trees[chrom].split_overlaps()
    for ovp in trees[chrom]:
        sys.stdout.write(chrom + "\t" + str(ovp[0]) + "\t" + str(ovp[1]) + "\n")
    
