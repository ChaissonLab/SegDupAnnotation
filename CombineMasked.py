#!/usr/bin/env python
import sys
fileList=sys.argv[1]
originalFai=sys.argv[2]
outputFile=sys.argv[3]

origFaiFile=open(originalFai)
origOrder = [ l.split()[0] for l in origFaiFile]


chroms={}

ifl=open(fileList)
inFiles=ifl.readlines()

for line in inFiles:


    tmp=open(line.rstrip())
    header=tmp.readline()[1:]
    vals=header.split("/")

    pre="/".join(vals[:-1])
    suf=int(vals[-1])
    tmp.close()
    
    if pre not in chroms:
        chroms[pre] = {}
    chroms[pre][suf] = line.rstrip()

outFile=open(outputFile,'w')


for origContig in origOrder:
    if origContig not in chroms:
        sys.stderr.write("ERROR: Missing contig: " + origContig + "\n")
        sys.exit(1)
    pre=origContig
    outFile.write(">" + pre + "\n")
    nParts=max(chroms[pre].keys())
    seq=""
    for i in range(0,nParts+1):
        partFile=open(chroms[pre][i])
        partFile.readline()
        seq+="".join([l.strip() for l in partFile.readlines()])
 
    print("chrom " + pre + " has length " + str(len(seq)))
    last=int(len(seq)/60)
    chromSeq="\n".join([seq[j*60:(j+1)*60] for j in range(0,last) ])
    if len(seq) % 60 != 0:
        chromSeq += "\n" + seq[last*60:]
    chromSeq += "\n"
    
    outFile.write(chromSeq)
        

