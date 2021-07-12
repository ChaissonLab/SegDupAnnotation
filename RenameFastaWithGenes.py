#!/usr/bin/env python
import sys
fastaFile=open(sys.argv[1])
bed12=open(sys.argv[2])

geneNames=[line.split()[3] for line in bed12]
index=0
for line in fastaFile:
    if line[0] == ">":
        sys.stdout.write(">" + geneNames[index] + "/" + line[1:].strip() + "\n")
        index+=1
    else:
        sys.stdout.write(line)
