#!/usr/bin/env python
import sys
from Bio import Seq
from Bio import SeqIO
from Bio import SeqRecord

inFile = open(sys.argv[1])
kmer = int(sys.argv[2])
index=0
for rec in SeqIO.parse(inFile, "fasta"):
    seq=str(rec.seq)
    for i in range(0,len(seq)-kmer+1):
        substr=seq[i:i+kmer]
        print(">"+str(index)  + "\n" + substr)

