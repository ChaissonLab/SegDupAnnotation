#!/usr/bin/env python
from Bio import Seq
from Bio import SeqIO
from Bio import SeqRecord
import sys

inFile = open(sys.argv[1])
L=int(sys.argv[2])
base=sys.argv[3]

for seqRec in SeqIO.parse(inFile, "fasta"):
    seq = str(seqRec.seq)
    seqLen = len(seq)

    nSeq = int(seqLen/L)
    if seqLen % L > 0:
        nSeq+=1
    for idx in range(0,nSeq):
        start=idx*L
        end=min((idx+1)*L, seqLen)
        sub = seq[start:end]
#        seqRec.id=seqRec.id.replace("/", "_").replace("|","_")
        outFile=open(base+ "." + seqRec.id + "_"+str(idx) + ".fasta", 'w')
        outFile.write(">"+seqRec.id+"/"+str(idx)+"\n")
        last=int(len(sub)/60)
        lines="\n".join([sub[j*60:(j+1)*60] for j in range(0,last)])
        
        if len(sub) % 60 > 0:
            lines += "\n" + sub[last*60:]
        outFile.write(lines)
        outFile.close()


