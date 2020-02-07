#!/usr/bin/env python
import sys
inFile = open(sys.argv[1])
for line in inFile:
    vals=line.split()
    cn=vals[3].split(",")
    i=0
    p=0
    while (i < len(cn)-1):
        if cn[i] != cn[i+1]:
            sys.stdout.write(vals[0] + "\t" + str(int(vals[1])+p*100) + "\t" + str(int(vals[1])+i*100) + "\t" + cn[i]+ "\n")
            p=i
        i+=1
    sys.stdout.write(vals[0] + "\t" + str(int(vals[1])+p*100) + "\t" + str(int(vals[1])+i*100) + "\t" + cn[i-1] + "\n")
            
