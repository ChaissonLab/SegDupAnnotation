#!/usr/bin/env python
import sys

for line in sys.stdin:
    if line[0] == "#":
        continue
    vals=line.split()

    sys.stdout.write("\t".join(vals[0:6]) + "\n")
    sys.stdout.write("\t".join(vals[3:6] + vals[0:3]) + "\n")


