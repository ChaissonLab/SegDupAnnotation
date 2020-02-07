#!/usr/bin/env python
import sys
import argparse
ap=argparse.ArgumentParser(description="Combine multiple rnaseq peaks into one file.")
ap.add_argument("--out", help="Output file", required=True)
ap.add_argument("--input", help="Input files", required=True, nargs="+")
args=ap.parse_args()

inFile = [open(v) for v in args.input ]
lines = [f.readlines() for f in inFile ]


outFile = open(args.out, 'w')



