#!/usr/bin/env python
import sys
fileList=sys.argv[1]
originalFai=sys.argv[2]
outputFile=sys.argv[3]
chromoutFilename=sys.argv[4]

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
chromOutFile=open(chromoutFilename,'w')
wroteChromOut=False
chromOutSeq=""
for origContig in origOrder:
    if origContig not in chroms:
        sys.stderr.write("ERROR: Missing contig: " + origContig + "\n")
        sys.exit(1)
    pre=origContig
    outFile.write(">" + pre + "\n")
    nParts=max(chroms[pre].keys())
    seq=""
    chromOutSeq=""
    for i in range(0,nParts+1):
        fn=chroms[pre][i]
        partFile=open(fn)
        
        outName=fn[:fn.rfind(".")] + ".out"
        cf=open(outName)
        outHeader = [cf.readline(), cf.readline(), cf.readline()]
        if (wroteChromOut is False):
            chromOutSeq="\n".join(outHeader)
            wroteChromOut = True

        for out in cf:
            vals=out.split()
            vals[5] = str(int(vals[5]) + len(seq))
            vals[6] = str(int(vals[6]) + len(seq))
            chromOutSeq+="\t".join(vals) + "\n"
            

        partFile.readline()
        seq+="".join([l.strip() for l in partFile.readlines()])
 
    print("chrom " + pre + " has length " + str(len(seq)))
    last=int(len(seq)/60)
    chromSeq="\n".join([seq[j*60:(j+1)*60] for j in range(0,last) ])
    if len(seq) % 60 != 0:
        chromSeq += "\n" + seq[last*60:]
    chromSeq += "\n"
    
    outFile.write(chromSeq)
    chromOutFile.write(chromOutSeq)
        

