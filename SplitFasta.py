#!/usr/bin/env python
from Bio import Seq
from Bio import SeqIO
from Bio import SeqRecord
import argparse
ap = argparse.ArgumentParser(description="Split a fasta file into small pieces")
ap.add_argument("fasta", help="Input fasta file")
ap.add_argument("L", help="Length of piece", type=int)
args=ap.parse_args()
inFile = open(args.fasta)
L=args.L

for seqRec in SeqIO.parse(inFile, "fasta"):
    seq = str(seqRec.seq)
    seqLen = len(seq)
    for i in range(0,int(seqLen/L)):
        start=i*L
        end=min((i+1)*L, seqLen)
        sub = seq[start:end]
        print(">"+seqRec.id + "/" +str(start) +"/"+ str(end) )
        print(sub)



