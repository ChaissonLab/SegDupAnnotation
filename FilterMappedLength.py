#!/usr/bin/env python
import pysam
import sys
af = pysam.AlignmentFile(sys.argv[1])
source=sys.argv[2]

strict=True
if source != "gencode":
    strict = False

outFile = sys.stdout

sys.stdout.write(str(af.header))
nFilt = 0
total = 0
for aln in af.fetch():
    astat=aln.get_cigar_stats()
    nHardClip=astat[0][5]
    nMatch = astat[0][0] + astat[0][7] + astat[0][8]
    if strict and "protein_coding" not in aln.query_name:
        continue
    if (aln.query_sequence is None):
        continue
    qLen = len(aln.query_sequence) + nHardClip

    if nMatch / qLen > 0.5:
        sys.stdout.write(aln.to_string() +"\n")

