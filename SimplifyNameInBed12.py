#!/usr/bin/env python
import sys
for line in sys.stdin:
    vals=line.split()
    if len(vals[3].split("|")) >= 9:
        vals[3] = vals[3].split("|")[5]
    else:
        vals[3] = vals[3].split("|")[0]
    sys.stdout.write("\t".join(vals) + "\n")
