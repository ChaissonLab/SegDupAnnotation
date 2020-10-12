#!/usr/bin/env python
import sys
import argparse
ap=argparse.ArgumentParser(description="Convert bed file of dups into circos files")
ap.add_argument("--bed", help="Input bed file", required=True)
ap.add_argument("--col", help="Bed file of genes in collapsed duplications. The name of the gene is in the 4th column")
ap.add_argument("--collapsed", help="Bed file of collaped duplications", required=True)
ap.add_argument("--links", help="Links file", default="resolved_dups.txt")
ap.add_argument("--labels", help="Labels for links", default="resolved_dups.labels.txt")

args=ap.parse_args()
linksf=open(args.links,'w')
labelsf=open(args.labels,'w')
bedf=open(args.bed)
colf=open(args.collapsed)

for line in bedf:
    vals=line.split()
    vals[12] = "var" + vals[12]
    sC=vals[12]
    sS=vals[13]
    sE=vals[14]
    vals[15] = "var" + vals[15]
    tC=vals[15]
    tS=vals[16]
    tE=vals[17]
    name=vals[3].split("|")[5]
    if sC!=tC:
        color="red"
#        color="set1-3-qual-1"
    else:
        color="blue"
#        color="set1-3-qual-2"
    line=vals[12:18] + ["color="+color]
    linksf.write("\t".join(line) + "\n")
    
    labelsf.write(sC + "\t" + sS + "\t" + sE + "\t" + name + "\n")

linksf.close()
labelsf.close()
            
