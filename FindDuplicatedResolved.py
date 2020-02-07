#!/usr/bin/env python
import sys


bedName = sys.argv[1]
dupName = sys.argv[2]

dupListFile = open(dupName)
dups = {l.strip() : True for l in dupListFile }
bedFile = open(bedName)
for line in bedFile:
    vals=line.split()
    if vals[3] in dups:
        sys.stdout.write(line)

