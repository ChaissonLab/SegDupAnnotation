#!/usr/bin/env python
import sys
refSeqFile=open("/home/cmb-16/nobackups/mjc/projects/VGP/rna-seq/Refseq2Gene.txt")
seqMap = { l.split()[0]:  l.split()[1] for l in refSeqFile }
for line in sys.stdin:
    vals=line.split()
    for i in range(0,len(vals)):
        if vals[i].split(".")[0] in seqMap:
            vals[i] = seqMap[vals[i].split(".")[0]]
    sys.stdout.write("\t".join(vals) + "\n")
