#!/usr/bin/env python
import sys
for line in sys.stdin:
    vals=line.split()
    vals[3] = vals[3].split("|")[5]
    sys.stdout.write("\t".join(vals) + "\n")
