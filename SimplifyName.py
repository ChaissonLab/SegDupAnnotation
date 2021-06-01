#!/usr/bin/env python
import sys
for line in sys.stdin:
    vals=line.split()
    for i in range(0,len(vals)):
        if vals[i].count("|") > 7:                   
            vals[i] = vals[i].split("|")[5]
    sys.stdout.write("\t".join(vals) + "\n")
