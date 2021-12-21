#!/usr/bin/env python
import sys
if len(sys.argv) < 3:
    print("Error: AppendOutFiles.py result.out first.out [second.out] ...")
    sys.exit(1)
    
outFile=open(sys.argv[1],'w')
f1=open(sys.argv[2])
lineNumber=0
maxIndex=0
for line in f1:
    outFile.write(line)
    lineNumber+=1
    if (lineNumber > 3):
        vals=line.split()
        index=vals[-1]
        if index == "*":
            index=vals[-2]
        index=int(index)
        if index > maxIndex:
            maxIndex=index
curIndex=maxIndex+1
for other in sys.argv[3:]:
    f=open(other)
    lineNumber=0
    for line in f:
        lineNumber +=1
        if lineNumber > 3:
            vals=line.split()
            index=vals[-1]
            if index == "*":
                index=vals[-2]
            index=int(index)
            index+=curIndex
            if (len(vals) < 14):
                print("ERROR with line")
                print(line)
                print(other)
                sys.exit(1)
            newLine="{:>5} {:>7.1} {:>5.1} {:>5.1}  {} {:>9} {:>8} {:>10} {:1} {:8} {:>23} {:>10} {:>7} {:>8} {:>10d}".format(vals[0],vals[1],vals[2],vals[3],
                                                                                                               vals[4], vals[5], vals[6], vals[7], vals[8],
                                                                                                               vals[9], vals[10], vals[11], vals[12],
                                                                                                               vals[13], index)
            if (vals[-1] == "*"):
                newLine+= " *"
            newLine += "\n"
            outFile.write(newLine)
