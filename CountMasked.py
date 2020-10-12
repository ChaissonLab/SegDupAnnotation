#!/usr/bin/env python
import sys
from Bio import Seq
from Bio import SeqIO
from Bio import SeqRecord
printCount=False
if len(sys.argv) > 2:
    if sys.argv[2] == "--count":
        printCount=True

for rec in SeqIO.parse(sys.argv[1], "fasta"):
    n = rec.seq.count("N") + rec.seq.count("a") + rec.seq.count("c") + rec.seq.count("g") + rec.seq.count("t") + rec.seq.count("n")
    output =rec.id + "\t{:2.2f}".format(n/len(rec.seq))
    if printCount:
        output+= "\t" + str(n)
    print(output)
