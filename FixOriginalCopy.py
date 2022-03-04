#!/usr/bin/env python
import sys
for line in sys.stdin:
    vals=line.strip().split()
    if "Original" in vals[-1]:
        vals[-1] = "Original"
    else:
        vals[-1] = "Copy"
    sys.stdout.write("\t".join(vals) + "\n")
