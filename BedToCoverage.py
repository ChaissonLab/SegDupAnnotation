#!/usr/bin/env python

import numpy as np
import argparse
import sys
ap = argparse.ArgumentParser(description="transform a bed into cvoerage.")

ap.add_argument("bed", help="Input bed file.")
ap.add_argument("bin", help="Bin size", type=int)
ap.add_argument("fai", help="Genome fai file.")
ap.add_argument("--out", help="Output file", default="/dev/stdout")
ap.add_argument("--minLength", help="Only consider alignments of this length", type=int, default=0)
ap.add_argument("--counts", help="Print bin counts, not whole bed.", default=False, action='store_true')
ap.add_argument("--totals", help="Print totals in bins, not coverage.", default=False, action='store_true')
args=ap.parse_args()

faiFile = open(args.fai)
bedFile = open(args.bed)
outFile = open(args.out, 'w')

fai = { line.split()[0]: int(line.split()[1]) for line in faiFile }

chroms = {}
for chrom in fai.keys():
    nBins = fai[chrom]/args.bin + np.ceil(float(fai[chrom]%args.bin)/args.bin) + 1
    if nBins > 0:
        chroms[chrom] = np.zeros(int(nBins), dtype=np.int32)

nReads = 0

for bedLine in bedFile:
    vals = bedLine.split()
    chrom = vals[0]
    if (chrom not in chroms):
        continue
    gStart = int(vals[1])
    gEnd   = int(vals[2])
    start = int(gStart/args.bin)
    end   = int(gEnd/args.bin + 1)
    if (end + 1 >= len(chroms[chrom])):
        continue    
    if (gEnd - gStart > args.minLength):
        if (args.totals == False):
            chroms[chrom][start] +=1
            chroms[chrom][end+1] -=1
    nReads +=1
ck = list(chroms.keys())
if (args.totals == False):

    for cki in ck:
        chrom = chroms[cki]
        cov = 0
        curInc = 0
        for i in range(0,len(chrom)-1):
            curInc = chrom[i]
            cov += curInc
            chrom[i] = cov

ck.sort()

if (args.counts):
    counts = {}
    for chrom in ck:
        chromBins = chroms[chrom]
        for i in range(0,len(chromBins)-1):
            if (chromBins[i] not in counts):
                counts[chromBins[i]] = 1
            else:
                counts[chromBins[i]]+=1
    for c,v in counts.iteritems():
        outFile.write("{}\t{}\n".format(c,v))
    sys.exit(0)


for chrom in ck:
    nBins = len(chroms[chrom])-1
    lines = ""
    chromBins = chroms[chrom]
    for i in range(0,nBins):
        lines += "{}\t{}\t{}\t{}\n".format(chrom, i*args.bin, min((i+1)*args.bin, fai[chrom]), chromBins[i])
        if ( (i+1)%100000 == 0 or i == nBins-1):
            outFile.write(lines)
            lines = ""
#            sys.stderr.write("{} {}\n".format(chrom, i+1))
