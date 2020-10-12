#!/usr/bin/env python
from Bio import SeqIO
from Bio import Seq
from Bio import SeqRecord
import sys
f=open(sys.argv[1])
for rec in SeqIO.parse(f,"fasta"):
    i=0
    s=str(rec.seq)
    sys.stderr.write("processing " + rec.id + "\n")
    while i < len(s):
        j=i+1
        while j < len(s) and s[j] != 'N' and s[j] != 'n':
            j+=1
        end=j


        while j < len(s) and (s[j] == 'N'or s[j] == 'n'):
            j+=1
        print(str(end-i))
        i=j+1


