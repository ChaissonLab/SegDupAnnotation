#!/usr/bin/env python

import sys
bedFile = open(sys.argv[1])
blacklistFile = open(sys.argv[2])

blacklist = {v.strip() : True  for v in blacklistFile.readlines()  }

for line in bedFile:
    vals=line.split()
    if vals[3] not in blacklist:
        sys.stdout.write(line)

