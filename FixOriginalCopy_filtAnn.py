#!/usr/bin/env python
import sys
for line in sys.stdin:
    vals=line.strip().split()
    if "Original" in vals[-2]:
        vals[-2] = "Original"
    else:
        vals[-2] = "Copy"
    sys.stdout.write("\t".join(vals) + "\n")
