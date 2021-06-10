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

    print(line)
    tmp=open(line.rstrip())
    try:
        header=tmp.readline()[1:]
    except ValueError:
        import pdb
        pdb.set_trace()

    h=header.split(" ")
    if len(h) != 2:
        import pdb
        pdb.set_trace()
    
    [chrom, annot]=h

        
    
    annotVals=annot.split("/")

    
    suf=int(annotVals[-1])
    tmp.close()
    
    if chrom not in chroms:
        chroms[chrom] = {}
    chroms[chrom][suf] = line.rstrip()
    
#print(chroms)

outFile=open(outputFile,'w')
chromOutFile=open(chromoutFilename,'w')
wroteChromOut=False
chromOutStr=""
repeatIndex=1
for origContig in origOrder:
    if origContig not in chroms:
        sys.stderr.write("ERROR: Missing contig: " + origContig + "\n")
        sys.exit(1)
    pre=origContig
    outFile.write(">" + pre + "\n")
    nParts=max(chroms[pre].keys())
    seq=""
    chromOutSeq=""

    allChromRepeats=[]
    toMask=[]
    for i in range(0,nParts+1):

#        print("pre " + pre)


        fn      =chroms[pre][i]
        partFile=open(fn)
        header=partFile.readline()
        headerVals = header.split()[1].split("/")
        seqStart   = int(headerVals[-4])
        seqEnd     = int(headerVals[-3])
        seqOffset  = int(headerVals[-2])
        seqIndex   = int(headerVals[-1])
        
        
        outName=fn[:fn.rfind(".")] + ".out"
        cf=open(outName)
        outHeader = [cf.readline(), cf.readline(), cf.readline()]
        if wroteChromOut is False:
            chromOutFile.write("".join(outHeader))
            wroteChromOut = True

#        import pdb
#        pdb.set_trace()
        for out in cf:
            vals=out.split()
            vals[5] = str(int(vals[5]) + seqStart)
            vals[6] = str(int(vals[6]) + seqStart)
            rep=[seqStart, seqEnd, seqOffset, seqIndex, int(vals[5]), int(vals[6]), vals]
            allChromRepeats.append(rep)

        maskedSeq="".join([l.strip() for l in partFile.readlines()])
        maskedSeq=maskedSeq[seqOffset:]

        seq+=maskedSeq

    # Now process all repeats to remove overlaps.
    i=0
    to_mask=[]
    while i < len(allChromRepeats):
        j=i
        # Set j to the first repeat of the overlapping sequence
        while j < len(allChromRepeats) and allChromRepeats[j][3] == allChromRepeats[i][3]:
            j+=1
        k=i
        c=j
        while k < j and (j == len(allChromRepeats) or allChromRepeats[k][4] < allChromRepeats[j][4]):
            to_mask.append((allChromRepeats[k][4], allChromRepeats[k][5]))
            chromOutFile.write("{:>5} {:>6} {:>4} {:>4}  {:20s}{:>9}{:>9}{:>11} {} {:<15} {:<18} {:>5} {:>6} {:>6} {:>5}\n".format(allChromRepeats[k][6][0], allChromRepeats[k][6][1], allChromRepeats[k][6][2], allChromRepeats[k][6][3], allChromRepeats[k][6][4], allChromRepeats[k][6][5],allChromRepeats[k][6][6],allChromRepeats[k][6][7], allChromRepeats[k][6][8],allChromRepeats[k][6][9],allChromRepeats[k][6][10],allChromRepeats[k][6][11],allChromRepeats[k][6][12],allChromRepeats[k][6][13],str(repeatIndex)))
            repeatIndex+=1
            k+=1
            
        print("For seq " + str(allChromRepeats[i][3]) + " copied until " + str(k) + " and next begins at " + str(j))
        i=j

    seqArray=list(seq)
    for i in range(0,len(to_mask)):
        for j in range(to_mask[i][0]-1, to_mask[i][1]-1):
            seqArray[j] = seqArray[j].lower()

    seq=''.join(seqArray)
    
    print("chrom " + pre + " has length " + str(len(seq)))
    last=int(len(seq)/60)
    chromSeq="\n".join([seq[j*60:(j+1)*60] for j in range(0,last) ])
    if len(seq) % 60 != 0:
        chromSeq += "\n" + seq[last*60:]
    chromSeq += "\n"
    
    outFile.write(chromSeq)
    chromOutFile.write(chromOutSeq)
        

