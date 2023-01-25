#!/usr/bin/env python
import sys
fastaFile=open(sys.argv[1])
bed12=open(sys.argv[2])
prefix=sys.argv[3]

geneNames=[(line.split()[3]).split("/")[0] for line in bed12]
index=0
for line in fastaFile:
    if line[0] == ">":
        sys.stdout.write(">" + geneNames[index] + ":" + prefix + "/" + line[1:].strip() + "\n")
        index+=1
    else:
        sys.stdout.write(line)

# Keon Rabbani
# Used to compare seg dups discovered in sda runs
